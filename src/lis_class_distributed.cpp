/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file lis_class.cpp
  \brief interface for linear solver built with for lis library.
*/

#include "lis_class_distributed.h"
#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdio>

/*
 
  Typical usage :

  linear_solver *solver = new lis ();
  solver->set_lhs_structure (lhs.num_owned_rows (), ir, jc);
  solver->analyze ();
  solver->set_lhs_data (xa);
  solver->factorize ();
  solver->set_rhs (rhs);
  if (solver->solver_type () == "iterative")
  solver->set_max_iterations (1000);
  solver->set_tolerance (1e-12);
  solver->solve ();
  solver->cleanup ();


*/

int
lis_distributed::init_lis_objects ()
{
  // Build lis structures
  row = new LIS_INT[n + 1];
  col = new LIS_INT[nnz];
  value = new LIS_SCALAR[nnz];

  for (int i = 0; i < nnz ; ++i)
    col[i] = jcol[i] - index_base;
   
  for (int i = 0; i < n + 1; ++i)
    row[i] = row_ptr[i] - row_ptr[0];

  if (initialized)
    {
      destroy_lis_objects ();
      initialized = false;
    }

  lis_matrix_create (MPI_COMM_WORLD, &A);
  lis_vector_create (MPI_COMM_WORLD, &b);
  lis_vector_create (MPI_COMM_WORLD, &x);
  lis_solver_create (&solver);
    
  lis_matrix_set_size (A, n, 0);  
  lis_matrix_set_csr (nnz, row, col, value, A);
  lis_vector_set_size (b, n, 0);
  lis_vector_duplicate (b, &x);
  LIS_INT ie;
  lis_vector_get_range (x, &row_s, &ie);
  assert (ie - row_s == n);
  
  initialized = true;
  return 1;
}

int
lis_distributed::assemble_lis_matrix ()
{
  for (int i = 0; i < nnz ; ++i)
    value[i] = data[i];
  
  lis_matrix_assemble (A);
  return 1;
}

int
lis_distributed::invoke_lis_solver ()
{

  for (int i = row_s; i < row_s + n; ++i)
    {
      lis_vector_set_value (LIS_INS_VALUE, i, rhs[i - row_s], b); 
      if (have_initial_guess)
        lis_vector_set_value (LIS_INS_VALUE, i, initial_guess[i - row_s], x);
    }

  if (! option_string_set)
    {
      std::stringstream opt;
      opt << "-maxiter " << max_iter
          << " -tol " << tolerance
          << " -i " << iterative_method
          << " -p " << preconditioner
          << " -conv_cond " << convergence_condition;

      if (have_initial_guess)
        opt << " -initx_zeros false ";
      else
        opt << " -initx_zeros true ";

      option_string = opt.str ();
      option_string_set = true;
    }

  auto options = std::unique_ptr<char> (new char[option_string.length () + 1]);
  std::copy (option_string.begin (),
             option_string.end (), options.get ());

  lis_solver_set_option (options.get (), solver);
  lis_solve (A, b, x, solver);

  lis_solver_get_iter (solver, &iter);
  lis_solver_get_time (solver, &time);

  //gather solution vector
  double temp = 0.0;
  for (int i = row_s; i < row_s + n; ++i)
    {
      lis_vector_get_value (x, i, &temp);
      rhs[i - row_s] = temp;
    }
  
  return 1;
}

void
lis_distributed::destroy_lis_objects ()
{  
  lis_solver_destroy (solver);
  lis_matrix_destroy (A);
  lis_vector_destroy (b);
  lis_vector_destroy (x);
}

void
lis_distributed::set_iterative_method (const std::string &s)
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
  option_string_set = false;
}


void
lis_distributed::set_preconditioner (const std::string &s)
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
  option_string_set = false;
}

void
lis_distributed::set_convergence_condition (const std::string &s)
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
  option_string_set = false;
}

void
lis_distributed::set_lhs_structure
(int n,
 std::vector<int> &ir,
 std::vector<int> &jc,
 matrix_format_t f)
{
  assert (f == csr);
  n_row = n;
  row_ptr.assign (n_row + 1, 0);
  jcol.assign (jc.size (), 0);

  row_ptr = ir;
  jcol = jc;
  nnz = jcol.size ();
}

int
lis_distributed::analyze () { return 1; }

void
lis_distributed::set_lhs_data (std::vector<double> &xa)
{ data = &*xa.begin (); }


int
lis_distributed::factorize ()
{   
  init_lis_objects ();
  int CHK = assemble_lis_matrix (); assert (CHK == 1);  
  return CHK;
}

int
lis_distributed::solve ()
{

  int retval = invoke_lis_solver ();

  if (verbose)
    std::cout << std::endl
              << "Number of iterations = " << iter
              << std::endl
              << "Elapsed time = " << time << std::endl;

  return retval;
}

