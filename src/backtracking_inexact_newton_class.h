/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file backtracking_inexact_newton_class.h
  \brief interface for a nonlinear solver.
*/

#ifndef HAVE_BACKTRACKING_INEXACT_NEWTON_H
#define HAVE_BACKTRACKING_INEXACT_NEWTON 1

#include <linear_solver.h>
#include <nonlinear_solver.h>
#include <abstract_nonlinear_problem.h>
#include <abstract_forcing_term.h>
#include <fstream>

class backtracking_inexact_newton : public nonlinear_solver
{
private :

  abstract_nonlinear_problem *problem;
  abstract_forcing_term *forcing;
  linear_solver *lin_solver;

  sparse_matrix lhs;
  std::vector<double> rhs;
  std::vector<double> *initial_guess;

  int max_iter;
  double min_residual;
  double tolerance;
  double forcing_value;

  double residual_norm;
  double step_norm;
  int iteration;

  norm_type norm_t;

  double t;
  double theta;
  double theta_min;
  double theta_max;  

  int verbose;
  std::string filename;
  std::ofstream fout;

  int rank, size;

public :

  /// Default costructor.
  /// The solution of each iteration will be print in output file
  /// if verbose_ = 2
  /// The solution of final iteration will be print in output file
  /// if verbose_ = 1
  /// Anyone solution will be print in output file
  /// if verbose_ = 0 (default)
  backtracking_inexact_newton (linear_solver *solver_, int verbose_ = 0) :
  nonlinear_solver ("Backtracking Inexact Newton"),
    verbose (verbose_)
  {
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    lin_solver = solver_;
    max_iter = 100;
    min_residual = 1e-10;
    tolerance = 1e-10;
    forcing_value = 1e-12;
    residual_norm = 0.0;
    step_norm = 0.0;
    iteration = 0;
    norm_t = L2;
    t = 10e-4;
    theta = 0;
    theta_min = 0;
    theta_max = 1;
    filename = "output.txt";
  };

  /// Set nonlinear problem to solve.
  void
  set_problem (abstract_nonlinear_problem *problem_)
  {
    problem = problem_;
  }

  /// Set forcing term.
  void
  set_forcing_term (abstract_forcing_term *forcing_)
  {
    forcing = forcing_;
  }

  /// Set initial guess.
  void
  set_initial_guess
  (std::vector<double> &initial_guess_)
  {
    initial_guess = &initial_guess_;
  }

  /// Solve system.
  int
  solve ();

  /// Set maximum number of nonlinear solver iterations.
  void
  set_max_iterations (int max_iter_)
  {
    max_iter = max_iter_;
  }

  /// Set tolerance of nonlinear solver.
  void
  set_tolerance (double tolerance_)
  {
    tolerance = tolerance_;
  }

  /// Set minimum residual of nonlinear solver.
  void
  set_min_residual (double min_residual_)
  {
    min_residual = min_residual_;
  }

  /// Set type of norm used by nonlinear solver.
  void
  set_norm_type (norm_type norm_t_)
  {
    norm_t = norm_t_;
  }

  /// Set backtracking parameters.
  /// \f$ t \in (0, 1) \f$
  /// \f$ \theta \in \[\theta_{min}, \theta_{max}\] \f$
  void
  set_backtracking_parameters
  (double t_,
   double theta_min_,
   double theta_max_)
  {
    t = t_;
    theta_min = theta_min_;
    theta_max = theta_max_;
  }
  /// \f$ \theta \f$ was chosen to minimize over
  /// \f$ \[\theta_{min}, \theta_{max} \f$ the quadratic
  /// \f$ p (\theta) \f$ for which \f$ p (0) = g (0), p' (0) = g' (0)\f$
  /// and \f$ p (1) = g (1)\f$, where
  /// \f$ g (\theta) = \|F (x_k + theta * s_k)\|_2^2 \f$
  /// \f$ a = \|F (x_k + s_k)\|_2^2 - \|F (x_k)\|^2_2 -
  /// 2*F (x_k)^TF' (x_k)s_k \f$
  /// \f$ b = 2F (x_k)^TF' (x_k)s_k \f$
  /// \f$ c = \|F (x_k)\|^2_2 
  void
  theta_choice (double a, double b, double c);

  /// Set maximum number of linear solver iterations.
  void
  set_max_iterations_of_linear_solver
  (int max_iteration)
  {
    lin_solver->set_max_iterations (max_iteration);
  }

  /// Set initial guess of linear solver.
  void
  set_initial_guess_of_linear_solver
  (std::vector<double> &initial_guess)
  {
    lin_solver->set_initial_guess (initial_guess);
  }

  /// Set initial tolerance of linear solver.
  void
  set_initial_tolerance_of_linear_solver
  (double initial_tolerance)
  {
    lin_solver->set_tolerance (initial_tolerance);
    forcing_value = initial_tolerance;
  }

  /// Set type of iterative method used by linear solver.
  void
  set_iterative_method_of_linear_solver
  (const std::string &iterative_method)
  {
    lin_solver->set_iterative_method (iterative_method);
  }

  /// Set type of preconditioner used by linear solver.
  void
  set_preconditioner_of_linear_solver
  (const std::string &preconditioner)
  {
    lin_solver->set_preconditioner (preconditioner);
  }

  /// Set type of convergence condition used by linear solver.
  void
  set_convergence_condition_of_linear_solver
  (const std::string &convergence_condition)
  {
    lin_solver->set_convergence_condition (convergence_condition);
  }

  /// Return name of linear solver used.
  const std::string&
  linear_solver_name ()
  {
    return lin_solver->solver_name ();
  }

  /// Return type of linear solver used.
  const std::string&
  linear_solver_type ()
  {
    return lin_solver->solver_type ();
  }

  /// Get the norm of the solution's residual
  void
  get_result_residual_norm (double &residual_norm_)
  {
    residual_norm_ =  residual_norm;
  }

  /// Get the solution found by nonlinear solver.
  void
  get_result_solution (std::vector<double> &solution)
  {
    solution = *initial_guess;
  }
  /// Get number of iterations when solve () ends
  void
  get_result_iterations (int &iterations_)
  {
    iterations_ = iteration;
  }

  /// Cleanup memory
  void
  cleanup ()
  {
    lin_solver->cleanup ();
  }

  /// Set the name of output file.
  void
  set_output_filename (const std::string &filename_)
  {
    filename = filename_;
  }
};

#endif
