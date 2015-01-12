/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <lis_operators.h>
#include <bim_sparse.h>
#include <stdlib.h>
#include <mpi.h>

void
lis_matrix_parallelization (sparse_matrix& sp,
                            sparse_matrix& sp_loc)
{
  int rank, size;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  int gn = 0;
  int n_loc = 0;
  if (rank == 0)
    gn = sp.size ();
  MPI_Bcast (&gn, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
    {
      sparse_matrix::col_iterator j;
      n_loc = gn / size + gn % size;
      sp_loc.resize (n_loc);

      for (int i = 0; i < n_loc; ++i)
        if (sp[i].size ())
          for (j = sp[i].begin (); j != sp[i].end (); ++j)
            sp_loc[i][sp.col_idx (j)] = sp.col_val (j);


      for (int k = 1; k < size; ++k)
        for (int i = n_loc + (k-1) * (gn / size);
             i < n_loc + k * (gn/size); ++i)
            {
              int row_size = sp[i].size ();
              MPI_Send (&row_size, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
              if (sp[i].size ())
                for (j = sp[i].begin (); j != sp[i].end (); ++j)
                  {
                    int col_idx = sp.col_idx (j);
                    double col_val = sp.col_val (j);
                    MPI_Send (&col_idx, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
                    MPI_Send (&col_val, 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
                  }
            }
    }
  else
    {

      MPI_Status *status;
      n_loc = gn / size;
      sp_loc.resize (n_loc);

      for (int i = 0; i < n_loc; ++i)
        {
          int row_size = 0;
          MPI_Recv (&row_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, status);
          if (row_size)
            {
              int col = 0;
              for (int j = 0; j < row_size; ++j)
                {
                  MPI_Recv (&col, 1, MPI_INT,
                            0, 0, MPI_COMM_WORLD, status);
                  MPI_Recv (&sp_loc[i][col], 1, MPI_DOUBLE,
                            0, 0, MPI_COMM_WORLD, status);
                }
            }
        }

    }
}

void
lis_vector_parallelization (std::vector<double>& v,
                            std::vector<double>& v_loc)
{
  int rank, size;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  int gn = 0;
  int n_loc = 0;
  if (rank == 0)
    gn = v.size ();
  MPI_Bcast (&gn, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank==0)
    {
      n_loc = gn / size + gn % size;
      v_loc.resize (n_loc);
      for (int i = 0; i < n_loc; ++i)
        v_loc[i] = v[i];

      for (int k = 1; k < size; ++k)
        for (int i = n_loc + (k - 1) * (gn / size);
             i < n_loc + k * (gn / size); ++i)
          MPI_Send (&v[i], 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Status *status;
      n_loc = gn / size;
      v_loc.resize (n_loc);
      for (int i = 0; i < n_loc; ++i)
        MPI_Recv (&v_loc[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, status);
    }
}

void
lis_vector_unification (LIS_VECTOR& v_loc,
                        std::vector<double>& v,
                        int vsize, int is, int ie)
{
  int rank, size;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  if (rank == 0)
    {
      v.resize (vsize);

      for (int i = is; i < ie; ++i)
        {
          double temp = 0.0;
          lis_vector_get_value (v_loc, i, &temp);
          v[i] = temp;
        }
      MPI_Status *status;
      for (int k = 1; k < size; ++k)
        {
          int is_loc, ie_loc;
          MPI_Recv (&is_loc, 1, MPI_INT, k, 0, MPI_COMM_WORLD, status);
          MPI_Recv (&ie_loc, 1, MPI_INT, k, 0, MPI_COMM_WORLD, status);
          double temp = 0.0;
          for (int i = is_loc; i < ie_loc; ++i)
            {
              MPI_Recv (&temp, 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, status);
              v[i] = temp;
            }
        }
    }
  else
    {
      MPI_Send (&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send (&ie, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      for (int i = is; i < ie; ++i)
        {
          double temp = 0.0;
          lis_vector_get_value (v_loc, i, &temp);
          MPI_Send (&temp, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
}

void
lis_solve_system (sparse_matrix& lhs,
                  std::vector<double>& rhs,
                  std::vector<double>& sol,
                  LIS_INT& iter, double& time,
                  int nnodes, linear_solver_option& option)
{
  int rank,size;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD ,&size);

  LIS_MATRIX A;
  LIS_VECTOR b;
  LIS_VECTOR x;

  LIS_SOLVER solver;

  std::vector<double> xa;
  std::vector<int>    ir, jc;
  int is, ie;

  lhs.csr (xa, jc, ir);

  int nnz = xa.size ();
  int n = lhs.size ();

  LIS_INT* row = (LIS_INT *)malloc ((n + 1) * sizeof (LIS_INT));
  LIS_INT* col = (LIS_INT *)malloc (nnz * sizeof (LIS_INT));
  LIS_SCALAR* value = (LIS_SCALAR *)malloc (nnz * sizeof (LIS_SCALAR));

  for (int i = 0;i < nnz; ++i)
    {
      col[i] = jc[i];
      value[i] = xa[i];
    }
  for (int i = 0; i < n + 1; ++i)
    {
      if (rank != 0)
        row[i] = ir[i] - ir[0];
      else
        row[i] = ir[i];
    }

  lis_matrix_create (LIS_COMM_WORLD, &A);
  lis_matrix_set_size (A, n, 0);
  lis_matrix_set_csr (nnz, row, col, value, A);

  lis_matrix_assemble (A);

  lis_vector_create (LIS_COMM_WORLD, &b);

  lis_vector_set_size (b, n, 0);
  lis_vector_get_range (b, &is, &ie);

  for (int i = is; i < ie; ++i)
    lis_vector_set_value (LIS_INS_VALUE, i, rhs[i - is], b);

  lis_vector_create (LIS_COMM_WORLD, &x);

  lis_vector_duplicate (b, &x);
  lis_solver_create (&solver);

  lis_solver_set_option (option.tolerance, solver);
  lis_solver_set_option (option.maxit,solver);
  lis_solver_set_option (option.other_opt,solver);

  lis_solve (A, b, x, solver);

  lis_solver_get_iter (solver, &iter);
  lis_solver_get_time (solver, &time);

  lis_vector_unification (x, sol, nnodes, is, ie);
  
  lis_solver_destroy (solver);
  lis_matrix_destroy (A);
  lis_vector_destroy (b);
  lis_vector_destroy (x);
}
