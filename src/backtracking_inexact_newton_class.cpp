/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file backtracking_inexact_newton_class.cpp
  \brief interface for a nonlinear solver.
*/

#include <abstract_nonlinear_problem.h>
#include <backtracking_inexact_newton_class.h>
#include <bim_sparse.h>
#include <operators.h>
#include <nonlinear_solver.h>
#include <linear_solver.h>

void
backtracking_inexact_newton::set_problem
  (abstract_nonlinear_problem *problem_)
{ problem = problem_; }

void
backtracking_inexact_newton::set_forcing_term
  (abstract_forcing_term *forcing_)
{ forcing = forcing_; }

void
backtracking_inexact_newton::set_initial_guess
  (std::vector<double> &initial_guess_)
{ initial_guess = &initial_guess_; }

int
backtracking_inexact_newton::solve ()
{

  std::vector<int> ir, jc;
  std::vector<double> xa;
  std::vector<double> lin_initial_guess;

  sparse_matrix mass_matrix;

  std::vector<double> f_old, f_new, df_gap;

if (rank == 0)
    {
      if (norm_t == L2 || norm_t == H1)
        {
          std::vector<double> ecoeff (problem->msh.nelements, 1.0);
          std::vector<double> ncoeff (problem->msh.nnodes, 1.0);
          bim3a_reaction (problem->msh, ecoeff, ncoeff, mass_matrix);
          if (norm_t == H1)
            bim3a_laplacian (problem->msh, ecoeff, mass_matrix);
        }

      problem->operator () (lhs, rhs, (*initial_guess));
      residual_norm = bim3a_norm2 (rhs);
    }

  MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (residual_norm < min_residual)
    {
      std::cout << "The Initial guess satisfies minimum residual."
                << std::endl;
      return 1;
    }

  if (rank == 0 && verbose >= 1)
    {
      fout.open (filename.c_str ());

      std::cout << "Result of Non Linear Test "
                << "\nwill be written in "
                << filename
                << std::endl << std::endl;
    }
  if (rank == 0)
    if (problem->msh.nnodes !=0)
      lin_initial_guess.assign (problem->msh.nnodes, 0.0);
    else
      lin_initial_guess.assign (1, 0.0);

  do
    {
      ++iteration;

      if (rank == 0 && verbose >= 1)
        std::cout << "\nNewton Iteration: "<< iteration << std::endl;

      if (rank == 0)
        {
          f_old.assign (rhs.size (), 0.0);
          for (unsigned int i = 0; i < rhs.size (); ++i)
            f_old[i] = - rhs[i];

          lhs.aij (xa, ir, jc, lin_solver->get_index_base ());
          lin_solver->set_lhs_structure (lhs.rows (), ir, jc);
        }

      lin_solver->analyze ();

      if (rank == 0)
        lin_solver->set_lhs_data (xa);

      lin_solver->factorize ();

      if (rank == 0)
        lin_solver->set_rhs (rhs);

      if (lin_solver->solver_type () == "iterative")
        {
          lin_solver->set_tolerance (forcing_value);
          lin_solver->set_initial_guess (lin_initial_guess);
        }

      lin_solver->solve ();

      if (lin_solver->solver_type () == "iterative")
        lin_initial_guess = rhs;

      if (rank == 0)
        {
          if (verbose == 2)
            fout << "Iteration: " << iteration << std::endl;

          std::vector<double> unew (initial_guess->size ());
          for (int i = 0; i < rhs.size (); ++i)
            unew[i] = rhs[i] + (*initial_guess)[i];

          double f_old_norm = residual_norm;
          double f_new_norm = 0.0;

          (*problem) (f_new, unew);
          f_new_norm = bim3a_norm2 (f_new);

          while (f_new_norm >
            (1 - t * (1 - forcing_value)) * f_old_norm)
            {
              std::vector<double> temp;
              double temp_norm = 0.0;

              bim3a_matrix_vector_product (lhs, rhs, temp);
              for (int i = 0; i < rhs.size (); ++i)
                temp_norm += f_old[i] * temp[i];

              double a = f_new_norm * f_new_norm -
                f_old_norm * f_old_norm - 2 * temp_norm;

              double b = 2* temp_norm;
              double c = f_old_norm * f_old_norm;

              theta_choice (a, b, c);

              for (unsigned int i = 0; i < rhs.size (); ++i)
                {
                  rhs[i] *= theta;
                  unew[i] = rhs[i] + (*initial_guess)[i];
                }

              forcing_value = 1 - theta * (1 - forcing_value);

              (*problem) (f_new, unew);

              f_new_norm = bim3a_norm2 (f_new);
            }

          bim3a_matrix_vector_product (lhs, rhs, df_gap);
          forcing_value = (*forcing)
            (f_old, f_new, df_gap, forcing_value);

          for (unsigned int i = 0; i < rhs.size (); ++i)
            (*initial_guess)[i] = rhs[i] + (*initial_guess)[i];

          if (norm_t == Inf)
            bim3a_norm (problem->msh, rhs, step_norm, norm_t);
          else
            {
              step_norm = 0.0;
              std::vector<double> temp;
              bim3a_matrix_vector_product (mass_matrix, rhs, temp);
              for (unsigned int i = 0; i < rhs.size (); ++i)
                step_norm += rhs[i] * temp[i];
              step_norm = sqrt (step_norm);
            }

          (*problem) (lhs, rhs, (*initial_guess));
          residual_norm = bim3a_norm2 (rhs);

         if (verbose == 2)
            for (unsigned int i = 0; i < initial_guess->size (); ++i)
              fout << (*initial_guess)[i] << std::endl;

          if (verbose >= 1)
            std::cout << "Step Error: " << step_norm << std::endl
                      << "Residual Error: " << residual_norm
                      << std::endl;
        }
      MPI_Bcast (&step_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&forcing_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  while (iteration <= max_iter
         && residual_norm > min_residual
         && step_norm > tolerance);

  if (rank == 0 && verbose == 1)
    {
      fout << "Solution: " << std::endl;
      for (unsigned int i = 0; i < initial_guess->size (); ++i)
        fout << (*initial_guess)[i] << std::endl;
    }

  if (rank == 0 && verbose >=1)
    fout.close ();

  if (iteration > max_iter)
    {
      if (rank == 0)
        std::cout << "Maximum number of iterations reached."
                  << std::endl
                  << "Nonlinear solver not converged."
                  << std::endl;
      return 0;
    }
  else
    {
      if (rank == 0)
        std::cout << "\nNonlinear solver converged."
                  << std::endl;
      return 1;
    }
}

void
backtracking_inexact_newton::theta_choice
(double a, double b, double c)
{
  double f_theta, f_min, f_max;
  theta = - b / (2 * a);
  f_theta = a * theta * theta + b * theta + c;
  f_min = a * theta_min * theta_min + b * theta_min + c;
  f_max = a * theta_max * theta_max + b * theta_max + c;
  if (theta > theta_max || theta < theta_min)
    if (f_min < f_max)
      theta = theta_min;
    else
      theta = theta_max;
  else
    if (f_theta <= f_min )
      if (f_theta >= f_max)
        theta = theta_max;
      else
        if (f_min <= f_max)
          theta = theta_min;
        else
          theta = theta_max;
}

void
backtracking_inexact_newton::set_max_iterations_of_linear_solver
  (int max_iteration)
{ lin_solver->set_max_iterations (max_iteration); }

void
backtracking_inexact_newton::set_initial_guess_of_linear_solver
  (std::vector<double> &initial_guess)
{ lin_solver->set_initial_guess (initial_guess); }

void
backtracking_inexact_newton::set_initial_tolerance_of_linear_solver
  (double initial_tolerance)
{
  lin_solver->set_tolerance (initial_tolerance);
  forcing_value = initial_tolerance;
}

void
backtracking_inexact_newton::set_iterative_method_of_linear_solver
  (const std::string &iterative_method)
{
  lin_solver->set_iterative_method (iterative_method);
}

void
backtracking_inexact_newton::set_preconditioner_of_linear_solver
  (const std::string &preconditioner)
{
  lin_solver->set_preconditioner (preconditioner);
}

void
backtracking_inexact_newton::set_convergence_condition_of_linear_solver
  (const std::string &convergence_condition)
{
  lin_solver->set_convergence_condition (convergence_condition);
}

void
backtracking_inexact_newton::get_result_residual_norm
  (double &residual_norm_)
{ residual_norm_ =  residual_norm; }

void
backtracking_inexact_newton::get_result_solution
  (std::vector<double> &solution)
{ solution = *initial_guess; }

void
backtracking_inexact_newton::get_result_iterations
  (int &iterations_)
{ iterations_ = iteration; }

void
backtracking_inexact_newton::cleanup ()
{ lin_solver->cleanup (); }

void
backtracking_inexact_newton::set_output_filename
  (const std::string &filename_)
{ filename = filename_; }
