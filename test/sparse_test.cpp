/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <bim_sparse.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>

int main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  if (rank == 0)
    {

      std::cout << "*****\ntest 1\n*****\n";

      sparse_matrix     sp;
      p_sparse_matrix   spp;

      std::vector<int>   i,j;
      std::vector<double>  a;

      std::vector<int>   rows,cols;

      rows.resize (2);
      rows[0] = 1; rows[1] = 2;

      cols.resize (3);
      cols[0] = 0; cols[1] = 1; cols[2] = 6;

      sp.resize (4);


      /*          col 0           col 1           col 2           col 3           col 4           col 5           col 6 */
      /* row 1 */ sp[1][0] = 1.0;                 sp[1][2] = 1.1; sp[1][3] = 1.3; sp[1][4] = 1.4;
      /* row 2 */ sp[2][0] = 2.0;                                                                                 sp[2][6] = 2.6;
      /* row 3 */                                                                                 sp[3][5] = 3.5;


      sp.extract_block_pointer (rows, cols, spp);

      spp.aij (a, i, j);
      std::cout << std::endl << "aij" << std::endl;

      for (unsigned int k = 0; k < sp.nnz; ++k)
        std::cout << i[k] << " "
                  << j[k] << " "
                  << a[k] << std::endl;

      std::cout << std::endl << "sparse" << std::endl;
      std::cout << spp;

    }

  MPI_Barrier (MPI_COMM_WORLD);

  sparse_matrix       lhs;
  std::vector<double> lrhs;
  std::vector<int>    ir, jc;
  std::vector<double> xa;
  std::vector<int>    rows2, cols2;
  p_sparse_matrix     llhs;

  if (rank == 0)
    {
      std::cout << "\n\n*****\ntest 2\n*****\n";

      lhs.resize (1000);
      for (int k = 0; k < 1000; ++k)
        lhs[k][k] = 2.0;

      for (int k = 0; k < 999; ++k)
        {
          lhs[k][k+1] = -1.0;
          lhs[k+1][k] = -1.0;
        }

      lrhs.resize (980);
      lrhs.assign (lrhs.size (), 1.0);

      rows2.resize (980);
      cols2.resize (980);

      for (int k = 0; k < 980; ++k)
        {
          rows2[k] = k+10;
          cols2[k] = k+10;
        }

      lhs.extract_block_pointer (rows2, cols2, llhs);
      llhs.aij (xa, ir, jc, 1);
    }

  mumps mumps_solver;

  if (rank == 0)
    mumps_solver.set_lhs_structure (llhs.rows (), ir, jc);
  mumps_solver.analyze ();

  // first solve
  if (rank == 0)
    mumps_solver.set_lhs_data (xa);

  mumps_solver.factorize ();

  if (rank == 0)
    mumps_solver.set_rhs (lrhs);

  mumps_solver.solve ();

  if (rank == 0)
    {
      std::cout << "\nresult of first solve\nwill be written in first_solve.txt\n";
      std::ofstream fout ("first_solve.txt");
      fout << std::endl;
      for (unsigned int k = 0; k < lrhs.size (); ++k)
        fout << lrhs[k] << std::endl;
      fout.close ();
    }

  // second solve
  if (rank == 0)
    {
      for (int k = 0; k < 1000; ++k)
        lhs[k][k] = 4.0;

      llhs.aij_update (xa, ir, jc, 1);

      mumps_solver.set_lhs_data (xa);
    }

  mumps_solver.factorize ();

  if (rank == 0)
    lrhs.assign (lrhs.size (), 1.0);

  if (rank == 0)
    mumps_solver.set_rhs (lrhs);

  mumps_solver.solve ();

  if (rank == 0)
    {
      std::cout << "\nresult of second solve\nwill be written in second_solve.txt\n";
      std::ofstream fout ("second_solve.txt");
      fout << std::endl;
      for (unsigned int k = 0; k < lrhs.size (); ++k)
        fout << lrhs[k] << std::endl;
      fout << std::endl;
      fout.close ();
    }

  mumps_solver.cleanup ();

  MPI_Finalize ();
  return (0);
}
