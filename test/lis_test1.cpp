/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <stdio.h>
#include <lis_config.h>
#include <lis.h>
#include <lis_operators.h>
#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <fstream>
#include <stdlib.h>

LIS_INT main (LIS_INT argc, char* argv[])
{
  int nnz, gn , n_loc, nnz_loc;
  int rank, size;
  sparse_matrix       sp, sp_loc;
  std::vector<int>    jc, ir;
  std::vector<double> xa;

  LIS_MATRIX A;
  LIS_INT *row,*col;
  LIS_SCALAR *value;

  lis_initialize (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  if (rank == 0)
    {
      std::cout << "\n\n*****\nLis Test 1\n*****\n";
      sp.resize (6);
      sp[0][0] = 10;
      sp[0][4] = -2;
      sp[1][0] = 3;
      sp[1][1] = 9;
      sp[1][5] = 3;
      sp[2][1] = 7;
      sp[2][2] = 8;
      sp[2][3] = 7;
      sp[3][0] = 3;
      sp[3][2] = 8;
      sp[3][3] = 7;
      sp[3][4] = 5;
      sp[4][1] = 8;
      sp[4][3] = 9;
      sp[4][4] = 9;
      sp[4][5] = 13;
      sp[5][1] = 4;
      sp[5][4] = 2;
      sp[5][5] = -1;

      sparse_matrix::col_iterator j;
      for (unsigned int i = 0; i < sp.size (); ++i){
        if (sp[i].size ())
          {
           for (j = sp[i].begin (); j != sp[i].end (); ++j)
             std::cout << sp.col_val (j) << " ";
           std::cout << std::endl;
          }
      }
      gn = sp.size ();
    }

  lis_matrix_parallelization (sp, sp_loc);

  sp_loc.csr (xa, jc, ir, 0);

  nnz_loc = xa.size ();
  n_loc = sp_loc.size ();

  row = (LIS_INT *) malloc ((n_loc + 1) * sizeof (LIS_INT));
  col = (LIS_INT *) malloc (nnz_loc * sizeof (LIS_INT));
  value = (LIS_SCALAR *) malloc (nnz_loc * sizeof (LIS_SCALAR));

  for (unsigned int i = 0; i < nnz_loc; ++i)
    {
      col[i] = jc[i];
      value[i] = xa[i];
    }
  for (int i = 0; i < n_loc + 1; ++i)
    {
      if (rank != 0)
        row[i] = ir[i] - ir[0];
      else
        row[i] = ir[i];
    }

  lis_matrix_create (LIS_COMM_WORLD, &A);
  lis_matrix_set_size (A, n_loc, 0);
  lis_matrix_set_csr (nnz_loc, row, col, value, A);

  LIS_INT error = lis_matrix_assemble (A);

  if (rank == 0)
    std::cout << "Matrix Assembled" << std::endl;
  if (rank == 0)
    {
      std::cout << "row: ";

      for (int i = 0; i < n_loc + 1; ++i)
        std::cout << row[i] << " ";

      std::cout << std::endl;
      std::cout << "col: ";

      for (int i = 0; i < nnz_loc; ++i)
        std::cout << col[i] << " ";

      std::cout << std::endl;
      std::cout << "val: ";

      for (int i = 0; i < nnz_loc; ++i)
        std::cout << value[i] << " ";

      std::cout << std::endl;

    }

  lis_matrix_destroy (A);
  lis_finalize ();
  return 0;
}
