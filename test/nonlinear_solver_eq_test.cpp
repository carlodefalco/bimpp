/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem:  
  \f[ x^2 - 2 = 0 \f]

  Exact Solution:
  \f[ x = \sqrt{2} \f]  
*/

#include <lis.h>
#include "mumps_class.h"
#include "lis_class.h"
#include "linear_solver.h"
#include "nonlinear_solver.h"
#include "operators.h"
#include "adaptive_inexact_newton_class.h"
#include "abstract_nonlinear_problem.h"
#include "abstract_forcing_term.h"
#include "forcing_class.h"
#include "bim_config.h"
#include "equation_class.h"

void
run_test_problem (nonlinear_solver *solver);

int main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  linear_solver *lis_solver = new lis ();
  //linear_solver *mumps_solver = new mumps ();

  nonlinear_solver *solver =
    new adaptive_inexact_newton (lis_solver, 2);

  run_test_problem (solver);

  MPI_Finalize ();
}

void
run_test_problem (nonlinear_solver *solver)
{
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  std::vector<double> exactsolution;

  std::vector<int>    bnodes;
  std::vector<double> vnodes;

  std::vector<double> uold;

  std::vector<double> fcoeff;

  double p = 2.0;

  abstract_nonlinear_problem *eq = new equation (p);
  abstract_forcing_term *forcing = new forcing_costant ();

  if (rank == 0)
    {
      std::cout << "\n\n*****\nNon Linear Equation Test\n*****\n";

      solver->set_output_filename ("NonLinear_Equation_Test.txt");

      std::cout << "Using solver named "
                << solver->solver_name ()
                << std::endl;

      exactsolution.resize (1);
      exactsolution[0] = sqrt (p);
      uold.resize (1);
      uold[0] = 0;

      eq->set_exact_solution (exactsolution);
    }

  if (rank == 0)
    {
      solver->set_problem (eq);
      solver->set_forcing_term (forcing);
      solver->set_initial_guess (uold);
    }

  solver->set_max_iterations (100);
  solver->set_tolerance (1e-10);
  solver->set_min_residual (1e-10);
  solver->set_norm_type (Inf);

  if (solver->linear_solver_type () == "iterative")
    {
      solver->set_max_iterations_of_linear_solver (1000);
      solver->set_iterative_method_of_linear_solver
              ("Conjugate Gradient");
      solver->set_initial_tolerance_of_linear_solver (0.5);
      solver->set_convergence_condition_of_linear_solver ("norm2_of_rhs");
    }

  bool converged = solver->solve ();

  if (rank == 0)
    {
      if (converged)
        {
          int iteration = 0;
          double residual_norm = 0.0;

          solver->get_result_iterations (iteration);
          solver->get_result_residual_norm (residual_norm);

      std::cout << "Solution: " << uold[0]
                    << std::endl
                    << "Total Newton's Iterations: "
                    << iteration << std::endl
                    << "Residual Norm: " << residual_norm << std::endl
                    << "Error: " << fabs (uold[0] - exactsolution[0])
                    << std::endl;
        }
      else
        {
          std::cerr << "Not Converged!" << std::endl;
          exit (-1);
        }
    }

  solver->cleanup ();
}
