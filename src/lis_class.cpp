/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file lis_class.cpp
  \brief interface for linear solver built with for lis library.
*/

#include <lis_class.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <cstring>

void
lis::set_lhs_structure
(int n,
 std::vector<int> &ir,
 std::vector<int> &jc,
 matrix_format_t f)
{
  n_row = n;
  row_ptr.clear ();
  row_ptr.resize (n_row + 1, 0);

  //aij_to_csr_format
  if (f == aij)
    {
      if (ir[0] == 1)
        for (int i = 0; i < ir.size (); ++i)
          {
            ir[i]--;
            jc[i]--;
          }
      for (int i = 0; i < ir.size (); ++i)
        row_ptr[ir[i]]++;

      for (int i = 0, cumsum = 0; i < n_row; ++i)
        {
          int temp = row_ptr[i];
          row_ptr[i] = cumsum;
          cumsum += temp;
        }
      row_ptr[n_row] = ir.size ();
    }
  else
    row_ptr = ir;

  jcol = jc;
}

int
lis::analyze ()
{
  //partitioning row_ptr and jcol
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  MPI_Bcast (&n_row, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
    {
      n = n_row / size;

      for (int k = 1; k < size; ++k)
        {
          i_s = row_ptr[n * k + n_row % size];
          nnz = row_ptr[n * (k + 1) + n_row % size] -
                row_ptr[n * k + n_row % size];

          MPI_Send (&i_s, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
          MPI_Send (&nnz, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
          MPI_Send (&jcol[i_s], nnz, MPI_INT, k, 0, MPI_COMM_WORLD);

          row_s = n * k + n_row % size;
          MPI_Send (&row_s, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
          MPI_Send (&row_ptr[row_s], n + 1,
                    MPI_INT, k, 0, MPI_COMM_WORLD);
        }
      n = n_row / size + n_row % size;
      nnz = row_ptr[n];
      i_s = 0;
      row_s = 0;
    }
  else
    {
      n = n_row / size;

      MPI_Status *status;
      MPI_Recv (&i_s, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, status);
      MPI_Recv (&nnz, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, status);

      jcol.resize (nnz, 0);
      row_ptr.resize (n + 1, 0);

      MPI_Recv (&jcol[0], nnz,
                MPI_INT, 0, 0, MPI_COMM_WORLD, status);
      MPI_Recv (&row_s, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, status);
      MPI_Recv (&row_ptr[0], n + 1,
                MPI_INT, 0, 0, MPI_COMM_WORLD, status);
    }
  return 1;
}

void
lis::set_lhs_data (std::vector<double> &xa)
{
  data = xa;
}

void
lis::set_rhs (std::vector<double> &rhs_)
{
  rhs = rhs_;
  rhs_it = rhs_.begin();
}

int
lis::solve ()
{
  //partitioning data and rhs
  int rank, size;
  LIS_INT iter;
  double time;
  LIS_SOLVER solver;
  LIS_MATRIX A;
  LIS_VECTOR b,x;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  if (rank == 0)
    {
      for (int k = 1; k < size; ++k)
        {
          n = n_row / size;
          nnz = row_ptr[n * (k + 1) + n_row % size] -
                row_ptr[n * k + n_row % size];
          i_s = row_ptr[n * k + n_row % size];
          row_s = n * k + n_row % size;

          MPI_Send (&data[i_s], nnz, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
          MPI_Send (&rhs[row_s], n, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
        }
      n = n_row / size + n_row % size;
      nnz = row_ptr[n];
      i_s = 0;
      row_s = 0;
    }
  else
    {
      MPI_Status *status;
      data.resize(nnz, 0);
      rhs.resize(n, 0);
      MPI_Recv (&data[0], nnz,
                MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, status);
      MPI_Recv (&rhs[0], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, status);
    }

  //Build lis structure and solve system
  LIS_INT *row = (LIS_INT *)malloc ((n + 1) * sizeof (LIS_INT));
  LIS_INT *col = (LIS_INT *)malloc (nnz * sizeof (LIS_INT));
  LIS_SCALAR *value = (LIS_SCALAR *)malloc (nnz * sizeof (LIS_SCALAR));

  for (int i = 0; i < nnz ; ++i)
    {
      col[i] = jcol[i];
      value[i] = data[i];
    }
  for (int i = 0; i < n + 1; ++i)
    row[i] = row_ptr[i] - row_ptr[0];

  lis_matrix_create (LIS_COMM_WORLD, &A);
  lis_matrix_set_size (A, n, 0);
  lis_matrix_set_csr (nnz, row, col, value, A);
  lis_matrix_assemble (A);

  lis_vector_create (LIS_COMM_WORLD, &b);
  lis_vector_set_size (b, n, 0);
  for(int i = row_s; i < row_s + n; ++i)
    lis_vector_set_value (LIS_INS_VALUE, i, rhs[i - row_s], b);

  lis_vector_create (LIS_COMM_WORLD, &x);
  lis_vector_duplicate (b, &x);

  lis_solver_create (&solver);

  std::stringstream opt;
  opt << "-maxiter " << max_iter
      << " -tol " << tolerance
      << " -i " << iterative_method
      << " -p " << preconditioner
      << " -conv_cond " << convergence_condition;

  std::string opt_ = opt.str ();
  char* options = new char[opt_.length () + 1];
  strcpy (options, opt_.c_str ());

  lis_solver_set_option (options, solver);

  lis_solve (A, b, x, solver);

  lis_solver_get_iter (solver, &iter);
  lis_solver_get_time (solver, &time);

  if (rank == 0)
    {
      std::cout << std::endl
                << "Number of iterations = " << iter
                << std::endl
                << "Elapsed time = " << time << std::endl;
    }

   delete[] options;

  //unificate solution vector
  if (rank == 0)
    {
      double temp = 0.0;
      for (int i = row_s; i < row_s + n; ++i)
        {
          lis_vector_get_value (x, i, &temp);
          *(rhs_it + i) = temp;
        }

      MPI_Status *status;
      for (int k = 1; k < size; ++k)
        {
          int loc_row_s, loc_n;
          MPI_Recv (&loc_row_s, 1,
                    MPI_INT, k, 0, MPI_COMM_WORLD, status);
          MPI_Recv (&loc_n, 1, MPI_INT, k, 0, MPI_COMM_WORLD, status);

          for (int i = loc_row_s; i < loc_row_s + loc_n; ++i)
            {
              MPI_Recv (&temp, 1, MPI_DOUBLE,
                        k, 0, MPI_COMM_WORLD, status);
              *(rhs_it + i) = temp;
            }
        }
    }
  else
    {
      double temp = 0.0;
      MPI_Send (&row_s, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send (&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      for (int i = row_s; i < row_s + n; ++i)
        {
          lis_vector_get_value (x, i, &temp);
          MPI_Send (&temp, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
  lis_solver_destroy (solver);
  lis_matrix_destroy (A);
  lis_vector_destroy (b);
  lis_vector_destroy (x);
  return 1;
}

void
lis::cleanup ()
{
  lis_finalize ();
}

void
lis::set_max_iterations (int max_iter_)
{
  max_iter = max_iter_;
}

void
lis::get_max_iterations (int &max_iter_)
{
  max_iter_ = max_iter;
}

void
lis::set_tolerance (double tol)
{
  tolerance = tol;
}

void
lis::get_tolerance (double &tol)
{
  tol = tolerance;
}

void
lis::set_iterative_method (const std::string &s)
{
  if (s == "conjugate_gradient")
    iterative_method = "cg";
  else if (s == "biconjugate_gradient")
    iterative_method = "bicg";
  else if (s == "bicg_stabilized")
    iterative_method = "bicgstab";
  else if (s == "jacobi")
    iterative_method = "jacobi";
  else if (s == "gauss_seidel")
    iterative_method = "gs";
  else if (s == "sor")
    iterative_method = "sor";
  else
    {
      std::cout << std::endl
                <<"Invalid Iterative Method"
                << std::endl
                << "Solve with default BiConjugate Gradient"
                << std::endl;
      iterative_method = "bicg";
    }
}

void
lis::get_iterative_method (std::string &s)
{
  s = iterative_method;
}

void
lis::set_preconditioner (const std::string &s)
{
  if (s == "none")
    preconditioner = "none";
  else if (s == "jacobi")
    preconditioner = "jacobi";
  else if (s == "ssor")
    preconditioner = "ssor";
  else
    {
      std::cout << std::endl
                <<"Invalid Preconditioner"
                << std::endl
                << "Solve without preconditioner"
                << std::endl;
      preconditioner = "none";
    }
}

void
lis::get_preconditioner (std::string &s)
{
  s = preconditioner;
}

void
lis::set_convergence_condition (const std::string &s)
{
  if (s == "norm2_of_residual")
    convergence_condition = "nrm2_r";
  else if (s == "norm2_of_rhs")
    convergence_condition = "nrm2_b";
  else
    {
      std::cout << std::endl
                <<"Invalid Convergence Condition"
                << std::endl
                << "Solve with default norm2_of_residual"
                << std::endl;
      convergence_condition = "nrm2_r";
    }
}

void
lis::get_convergence_condition (std::string &s)
{
  s = convergence_condition;
}
