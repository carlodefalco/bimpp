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

int
backtracking_inexact_newton::solve ()
{

  std::vector<int> ir, jc;
  std::vector<double> xa;
  std::vector<double> lin_initial_guess;

  if (rank == 0)
    {
      problem->operator () (lhs, rhs, (*initial_guess));
      bim3a_norm (problem->msh, rhs, residual_norm, norm_t);
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
    lin_initial_guess.assign (problem->msh.nnodes, 0.0);

  do
    {
      ++iteration;

      if (rank == 0 && verbose == 2)
        std::cout << "\nNewton Iteration: "<< iteration << std::endl;

      if (rank == 0)
        {
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

          std::vector<double> f_new, f_old;
          double f_old_norm = 0.0;
          double f_new_norm = 0.0;
          (*problem) (f_old, (*initial_guess));
          (*problem) (f_new, unew);
          bim3a_norm (problem->msh, f_old, f_old_norm, norm_t);
          bim3a_norm (problem->msh, f_new, f_new_norm, norm_t);

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
              bim3a_norm (problem->msh, f_new, f_new_norm, norm_t);
            }

          forcing_value = (*forcing) (problem, (*initial_guess),
                                      rhs, forcing_value, norm_t);

          for (unsigned int i = 0; i < rhs.size (); ++i)
            (*initial_guess)[i] = rhs[i] + (*initial_guess)[i];

          bim3a_norm (problem->msh, rhs, step_norm, norm_t);

          (*problem) (lhs, rhs, (*initial_guess));
          bim3a_norm (problem->msh, rhs, residual_norm, norm_t);
          if (verbose == 2)
            {
              for (unsigned int i = 0; i < initial_guess->size (); ++i)
                fout << (*initial_guess)[i] << std::endl;

              std::cout << "Step Error: " << step_norm << std::endl
                        << "Residual Error: " << residual_norm
                        << std::endl;
            }
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

      std::cout << "Residual Error: " << residual_norm
                << std::endl;
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
