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
  row_ptr.assign (n_row + 1, 0);

  jcol.assign (jc.size (), 0);

  //aij_to_csr_format
  if (f == aij)
    {
      ordering_map.assign (jc.size (), 0);
      for (unsigned int i = 0; i < ir.size (); ++i)
        row_ptr[ir[i] - index_base]++;

      for (unsigned int i = 0, cumsum = index_base; i < n_row; ++i)
        {
          int temp = row_ptr[i];
          row_ptr[i] = cumsum;
          cumsum += temp;
        }
      row_ptr[n_row] = ir.size () + index_base;
      for (unsigned int i = 0; i < jc.size (); ++i)
        {
          int row = ir[i];
          int dest = row_ptr[row];

          jcol[dest] = jc[i];
          ordering_map[i] = dest;

          row_ptr[row]++;
        }
      for (unsigned int i = 0, last = index_base; i <= n_row; ++i)
        {
          int temp = row_ptr[i];
          row_ptr[i] = last;
          last = temp;
        }
    }
  else
    {
      row_ptr = ir;
      jcol = jc;
    }
}

int
lis::analyze ()
{

  //partitioning row_ptr and jcol

  MPI_Bcast (&n_row, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
    {
      map_i_s.assign (size, 0);
      map_row_s.assign (size, 0);
      map_n.assign (size, 0);
      map_nnz.assign (size, 0);

      map_n[0] = n_row / size + n_row % size;
      map_nnz[0] = row_ptr[map_n[0]] - index_base;

      for (unsigned int k = 1; k < size; ++k)
        {
          map_n[k] = n_row / size;
          map_i_s[k] = row_ptr[map_n[k] * k + n_row % size] -
                       index_base;
          map_nnz[k] = row_ptr[map_n[k] * (k + 1) + n_row % size] -
                       row_ptr[map_n[k] * k + n_row % size];
          map_row_s[k] = map_n[k] * k + n_row % size;

        }
    }

  MPI_Scatter (&map_i_s[0], 1, MPI_INT,
               &i_s, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter (&map_row_s[0], 1, MPI_INT,
               &row_s, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter (&map_n[0], 1, MPI_INT,
               &n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter (&map_nnz[0], 1, MPI_INT,
               &nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
    MPI_Scatterv (&jcol[0], &map_nnz[0], &map_i_s[0], MPI_INT,
                  MPI_IN_PLACE, 0, MPI_INT, 0, MPI_COMM_WORLD);
  else
    {
      jcol.assign (nnz, 0);
      MPI_Scatterv (&jcol[0], &map_nnz[0], &map_i_s[0], MPI_INT,
        &jcol[0], nnz, MPI_INT, 0, MPI_COMM_WORLD);
    }
  //(commento da rimuovere in seguito)
  //non posso utilizzare scatter per row_ptr perchè ogni processo
  //riceve anche il primo valore che viene inviato al processo successivo
  if (rank == 0)
    for (unsigned int k = 1; k < size; ++k)
      MPI_Send (&row_ptr[map_row_s[k]], map_n[k] + 1, MPI_INT,
        k, 0, MPI_COMM_WORLD);
  else
    {
      MPI_Status *status = NULL;
      row_ptr.assign (n + 1, 0);
      MPI_Recv (&row_ptr[0], n + 1, MPI_INT,
        0, 0, MPI_COMM_WORLD, status);
    }
  return 1;
}

void
lis::set_lhs_data (std::vector<double> &xa)
{
  //(commento da rimuovere in seguito)
  //per non modificare xa penso sia necessario allocare nuova memoria
  if (ordering_map.size () != 0)
    {
      data = new double [xa.size ()];
      for (unsigned int i = 0; i < xa.size (); ++i)
      data[ordering_map[i]] = xa[i];
    }
  else
    data = &*xa.begin ();
}

void
lis::set_rhs (std::vector<double> &rhs_)
{
  rhs = &*rhs_.begin ();
}

void
lis::set_initial_guess (std::vector<double> &initial_guess_)
{
  have_initial_guess = true;
  initial_guess = &*initial_guess_.begin ();
}

int
lis::solve ()
{
  LIS_INT iter;
  double time;
  LIS_SOLVER solver;
  LIS_MATRIX A;
  LIS_VECTOR b,x;

  //partitioning data, rhs and initial_guess if exist
  MPI_Bcast (&have_initial_guess, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
    {
      MPI_Scatterv (&data[0], &map_nnz[0], &map_i_s[0], MPI_DOUBLE,
        MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Scatterv (&rhs[0], &map_n[0], &map_row_s[0], MPI_DOUBLE,
        MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (have_initial_guess)
        MPI_Scatterv (&initial_guess[0], &map_n[0], &map_row_s[0],
          MPI_DOUBLE, MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  else
    {
      data = new double[nnz];
      rhs = new double[n];

      MPI_Scatterv (&data[0], &map_nnz[0], &map_i_s[0], MPI_DOUBLE,
        &data[0], nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Scatterv (&rhs[0], &map_n[0], &map_row_s[0], MPI_DOUBLE,
        &rhs[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (have_initial_guess)
        {
          initial_guess = new double[n];
          MPI_Scatterv (&initial_guess[0], &map_n[0], &map_row_s[0],
            MPI_DOUBLE, &initial_guess[0], n, MPI_DOUBLE,
            0, MPI_COMM_WORLD);
        }
    }

  //Build lis structure and solve system
  LIS_INT *row = (LIS_INT *)malloc ((n + 1) * sizeof (LIS_INT));
  LIS_INT *col = (LIS_INT *)malloc (nnz * sizeof (LIS_INT));
  LIS_SCALAR *value = (LIS_SCALAR *)malloc (nnz * sizeof (LIS_SCALAR));

  for (unsigned int i = 0; i < nnz ; ++i)
    {
      col[i] = jcol[i] - index_base;
      value[i] = data[i];
    }
  for (unsigned int i = 0; i < n + 1; ++i)
    row[i] = row_ptr[i] - row_ptr[0];

  lis_matrix_create (LIS_COMM_WORLD, &A);
  lis_matrix_set_size (A, n, 0);
  lis_matrix_set_csr (nnz, row, col, value, A);
  lis_matrix_assemble (A);

  lis_vector_create (LIS_COMM_WORLD, &b);
  lis_vector_set_size (b, n, 0);
  lis_vector_create (LIS_COMM_WORLD, &x);
  lis_vector_duplicate (b, &x);
  for (unsigned int i = row_s; i < row_s + n; ++i)
    {
      lis_vector_set_value (LIS_INS_VALUE, i, rhs[i - row_s], b);
      if (have_initial_guess)
        lis_vector_set_value (LIS_INS_VALUE, i,
                              initial_guess[i - row_s], x);
    }

  lis_solver_create (&solver);

  std::stringstream opt;
  opt << "-maxiter " << max_iter
      << " -tol " << tolerance
      << " -i " << iterative_method
      << " -p " << preconditioner
      << " -conv_cond " << convergence_condition;
  if (have_initial_guess)
    opt << " -initx_zeros false ";
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

  delete [] options;

  //unificate solution vector
  double temp = 0.0;
  for (unsigned int i = row_s; i < row_s + n; ++i)
    {
      lis_vector_get_value (x, i, &temp);
      rhs[i - row_s] = temp;
    }
  if (rank == 0)
    MPI_Gatherv (MPI_IN_PLACE, 0, MPI_DOUBLE, &rhs[0],
      &map_n[0], &map_row_s[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv (&rhs[0], n, MPI_DOUBLE, &rhs[0],
      &map_n[0], &map_row_s[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  lis_solver_destroy (solver);
  lis_matrix_destroy (A);
  lis_vector_destroy (b);
  lis_vector_destroy (x);
  return 1;
}

void
lis::cleanup ()
{
  if (rank != 0)
    {
      delete [] data;
      delete [] rhs;
      if (have_initial_guess)
        delete [] initial_guess;
    }
  else if (ordering_map.size () != 0)
    delete [] data;

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
  if (s == "Conjugate Gradient")
    iterative_method = "cg";
  else if (s == "Biconjugate Gradient")
    iterative_method = "bicg";
  else if (s == "Bicg Stabilized")
    iterative_method = "bicgstab";
  else if (s == "Jacobi")
    iterative_method = "jacobi";
  else if (s == "Gauss Seidel")
    iterative_method = "gs";
  else if (s == "SOR")
    iterative_method = "sor";
  else
    {
      if (rank == 0)
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
  else if (s == "ilu")
    preconditioner = "ilu";
  else if (s == "ilut")
    preconditioner = "ilut";
  else if (s == "iluc")
    preconditioner = "iluc";
  else
    {
      if (rank == 0)
        std::cout << std::endl
                  <<"Invalid Preconditioner"
                  << std::endl
                  << "Solve with ilu[0]"
                  << std::endl;
      preconditioner = "ilu";
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
      if (rank == 0)
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
