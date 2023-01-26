/*
  Copyright (C) 2023 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/


#include <bim_sparse_distributed.h>
#include <bim_distributed_vector.h>
#include <lis_class_distributed.h>

//!     Assemble a discrete Jacobian on the mesh depicted
//!      below. 'X' marks vertices and elements owned by
//!      rank 0, while 'O' marks vertices and elements
//!      owned by rank 1. Ordering is column-major with the
//!      numbering shown below.
//!
//!      4       X---X---X---X---X---O---O---O---O---O---O
//!              |   |   |   |   |   |   |   |   |   |   |
//!              | X | X | X | X | X | O | O | O | O | O |
//!      3       X---X---X---X---X---O---O---O---O---O---O
//!              |   |   |   |   |   |   |   |   |   |   |
//!              | X | X | X | X | X | O | O | O | O | O |
//!      2       X---X---X---X---X---X---O---O---O---O---O
//!              |   |   |   |   |   |   |   |   |   |   |
//!              | X | X | X | X | X | O | O | O | O | O |
//!      1       X---X---X---X---X---X---O---O---O---O---O
//!              |   |   |   |   |   |   |   |   |   |   |
//!              | X | X | X | X | X | O | O | O | O | O |
//!      0       X---X---X---X---X---X---O---O---O---O---O
//!
//!              0   1   2   3   4   5   6   7   8   9   10

int
conn (int inode, int iel)
{
  int res = 0;
  int elcol, elrow;
  elcol = iel / 4;
  elrow = iel - elcol * 4;
  res = elcol * 5 + elrow;
  if (inode == 1) res += 1;
  if (inode == 2) res += 5;
  if (inode == 3) res += 6;
  return (res);
};

int
main (int argc, char *argv[])
{
  MPI_Init (&argc, &argv);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  int is, ie, is_elems, ie_elems, rank_owned;
  distributed_sparse_matrix A;


  if (rank == 0)
    {
      is_elems = 0;
      ie_elems = 20;
      is = 0;
      ie = 28;
      rank_owned = ie - is;
    }
  else if (rank == 1)
    {
      is_elems = 20;
      ie_elems = 40;
      is = 28;
      ie = 55;
      rank_owned = ie - is;
    }
  else
    {
      is_elems = 40;
      ie_elems = 40;
      is = 55;
      ie = 55;
      rank_owned = ie - is;
    }

  if (size == 1)
    {
      std::cout << "running serially\n" ;
      is_elems = 0;  is = 0;
      ie_elems = 40; ie = 55;
      rank_owned = ie - is;
    }


  //A.set_ranges (is, ie);
  A.set_ranges (rank_owned);
  distributed_vector b(rank_owned);
  distributed_vector x0(rank_owned);

  b.get_owned_data ().assign (rank_owned, 1.0);
  x0.get_owned_data ().assign (rank_owned, rank);

  std::cout << "rank " << rank << " size " << A.size ()
	    << " is " << is << " ie " << ie
	    << " is_elems " << is_elems
	    << "  ie_elems " << ie_elems << std::endl;

  std::vector<std::vector<double>> locmatrix =
    {{2,-1,-1,0}, {-1,2,0,-1},{-1,0,2,-1},{0,-1,-1,2}};


  double tmp;
  for (int iel = is_elems; iel < ie_elems; ++iel)
    for (int inode = 0; inode < 4; ++inode)
      for (int jnode = 0; jnode < 4; ++jnode)
	{
	  tmp = locmatrix[inode][jnode];
	  if (tmp != 0.)
	    A[conn (inode,iel)][conn (jnode,iel)] +=
	      tmp;
	}

  if (rank == 0)
    {
      A[0][0] += 1000; b(0) += 1000;
      A[54][54] += 1000; b(54) += 1000;
    }
  A.assemble ();
  b.assemble ();
  x0.assemble ();

  for (int ii = 0; ii < size; ++ii)
    {
      if (ii == rank)
	{
	  std::cout << "## rank " << rank << std::endl;
	  std::cout << A << std::endl;
	  std::cout << "\n\n";
	}
      MPI_Barrier (MPI_COMM_WORLD);
    }

  lis_distributed lis_solver;
  lis_solver.set_convergence_condition (std::string ("norm2_of_rhs"));
  lis_solver.set_tolerance (1.0e-17);

  std::vector<double> vals;
  std::vector<int> irow, jcol;

  A.csr (vals, jcol, irow, lis_solver.get_index_base ());

  // for (int ii = 0; ii < size; ++ii) {
  //   if (ii == rank) {
  //     std::cout << "rank = " << rank
  //		<< "rank_owned = " << rank_owned
  //		<< "vals.size() = " << vals.size()
  //		<< "jcol.size() = " << jcol.size()
  //		<< "irow.size() = " << irow.size()
  //		<< std::endl;
  //   }
  //   MPI_Barrier(MPI_COMM_WORLD);
  // }

  lis_solver.set_lhs_structure (A.rows (), irow, jcol);
  lis_solver.set_lhs_data (vals);

  lis_solver.set_rhs (b.get_owned_data ());
  lis_solver.set_initial_guess (x0.get_owned_data ());

  // Solve.
  lis_solver.analyze ();
  lis_solver.factorize ();
  lis_solver.solve ();

  if (rank == 0)
    std::cout << "x = [" << "\n";

  for (int jj = 0; jj < size; ++jj) {
    if (rank == jj)
      for (std::size_t ii = b.get_range_start (); ii < b.get_range_end (); ++ii) {
	std::cout << b[ii] << "\n";
      }
    MPI_Barrier (MPI_COMM_WORLD);
  }

  if (rank == 0)
    std::cout << "]:" << "\n";

  lis_solver.cleanup ();

  MPI_Finalize ();
  return 0;
}
