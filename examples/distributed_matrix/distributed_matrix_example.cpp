/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/


#include <bim_sparse_distributed.h>


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
  
  A.resize (55);


  if (rank == 0)
    {
      is_elems = 0;
      ie_elems = 20;
      is = 0;
      ie = 28;
      rank_owned = ie - is;
      std::cout << "if (rank == 0 || size == 1)" << std::endl;
    }
  else if (rank == 1)
    {
      is_elems = 20;
      ie_elems = 40;
      is = 28;
      ie = 55;
      rank_owned = ie - is;
      std::cout << "if (rank == 1 || size == 1)" << std::endl;
    }
  else
    {
      is_elems = 40;
      ie_elems = 40;
      is = 55;
      ie = 55;
      rank_owned = ie - is;
      std::cout << "else" << std::endl;
    }

  if (size == 1)
    {
      std::cout << "runing serially\n" ;
      is_elems = 0;  is = 0;
      ie_elems = 40; ie = 55;
      rank_owned = ie - is;
    }
  
  // A.set_ranges (is, ie);
  A.set_ranges (rank_owned);
  std::cout << "rank " << rank << " is " << is << " ie " << ie << " is_elems " << is_elems << "  ie_elems " << ie_elems << std::endl;  
        
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

  A.assemble ();
  
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
   
  MPI_Finalize ();
  return 0;
}
