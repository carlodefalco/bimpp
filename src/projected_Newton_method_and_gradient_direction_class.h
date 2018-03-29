/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file backtracking_inexact_newton_example8_class.h
  \brief interface for a nonlinear solver based on projected Newton method combined with prejected gradiet direction.
*/

#ifndef HAVE_PROJECTED_NEWTON_METHOD_AND_GRADIENT_DIRECTION_H
#define HAVE_PROJECTED_NEWTON_METHOD_AND_GRADIENT_DIRECTION_H 1


#include <linear_solver.h>
#include <nonlinear_solver.h>
#include <abstract_nonlinear_problem.h>
#include <abstract_forcing_term.h>
#include <fstream>

/// \brief Specific interface's class for a nonlinear solver.
class projected_Newton_method_and_gradient_direction : public nonlinear_solver
{
private :
  /// Pointer to the nonlinear problem used by nonlinear solver.
  abstract_nonlinear_problem *problem;

  /// Pointer to the forcing term used by nonlinear solver.
  abstract_forcing_term *forcing;

  /// Left Hand Side of the nonlinear problem linearized.
  sparse_matrix lhs;

  /// Right Hand Side of the nonlinear problem linearize.
  std::vector<double> rhs;

  /// Pointer to the initial guess of nonlinear solver.
  std::vector<double> *initial_guess;

  /// \brief Maximum number iterations of nonlinear solver.
  /// \details The iteration of nonlinear solver stops
  /// if iteration > max_iter [default = 100].
  int max_iter;

  /// \brief Minimum residual norm.
  /// \details Nonlinear solver stops if
  /// \f$ ||F(x)|| < min\_residual \f$ [default = 1e-10].
  double min_residual;

  /// \brief Tolerance for two successive iterations.
  /// \details The iteration of nonlinear solver stops
  /// if \f$||x_{new} - x_{old}|| < tolerance \f$
  /// [default = 1e-10].
  double tolerance;

  /// \brief Forcing value of linear solver.
  /// \details The linear solver iteration stops if
  /// convergence condition is satisfies.
  /// [default = 0.765518617913987]
  double forcing_value;

  /// \brief The norm of the residual in a specific nonlinear iteration.
  /// \details At the end of method solve ()
  /// it's the norm of the solution's residual.
  double residual_norm;

  /// \brief The norm of difference between two nonlinear iteration.
  double step_norm;

  /// \brief The iteration of nonlinear solver.
  /// \details At the and of method solve ()
  /// it's the number of iterations used by nonlinear solver.
  int iteration;

  /// \brief The type of norm used by nonlinear solver.
  /// \details [default = L2].
  norm_type norm_t;

  /// \brief Backtracking parameter.
  /// \details \f$||F (x_{new})|| > [1 - t (1-\eta_k)]||F (x_{old})||\f$
  /// [default = 1e-4].
  double t;

  /// \brief Backtracking parameter.
  /// \details \f$\Theta(\mathcal{P}(x_{new})) \leq \Theta(x_{old}) + 
  /// \sigma \nabla \Theta(x_{old})^T(\mathcal{P}(x_{new})-x_{old})\f$,
  /// with \f$ \Theta(x) = \frac{1}{2} ||F(x)||\f$
  /// [default = 1e-4].
  double sigma; 

  /// \brief Backtracking parameter for PN.
  /// \details Backtracking quantity
  /// computed for projected Newton direction.
  /// [default = 0.5].
  double theta;
  

  /// \brief Backtracking parameter for PG.
  /// \details Backtracking quantity
  /// computed for the projected gradient direction.
  /// [default = 0.8].
  double thetaPG;

  /// \brief Backtracking parameter.
  /// \details Lower bound of \f$ \theta \f$ [default = 0].
  double theta_min;

  /// \brief Backtracking parameter.
  /// \details Upper bound of \f$ \theta \f$ [default = 1].
  double theta_max;

  /// \brief Backtracking max iterations.
  /// \details maximum number of iterations in the backtracking while.
  /// [default = 20].
  double max_back_it;


  /// \brief The type of verbose.
  /// \details The solution of each iteration
  /// will be print in output file if verbose = 2;
  ///
  /// The solution of final iteration
  /// will be print in output file
  /// if verbose_ = 1 [default];
  ///
  /// Anyone solution will be print in output file
  /// if verbose_ = 0.
  
  int verbose;
  
  /// \brief Bounds for the solution.
  /// \details b1 lower bound, b2 upper bound.
  std::vector<double> b1;
  std::vector<double> b2;
  
  std::string filename;
  std::ofstream fout;

  int rank, size;

public :

  /// Pointer to the linear solver used by nonlinear solver.
  linear_solver *lin_solver;

  /// Default costructor.
  projected_Newton_method_and_gradient_direction
  (linear_solver *solver_, std::vector<double> &b1_ ,
   std::vector<double> &b2_ ,int verbose_ = 2, double t_ = 1e-4,
  double sigma_ = 1e-4, double theta_min_ =0, double theta_max_ =1, 
  double theta_ =0.5, double thetaPG_ =0.8, int max_back_it_=20) :
    nonlinear_solver ("projected_Newton_method_and_gradient_direction"),
    lin_solver (solver_),
    max_iter (100),
    min_residual (1e-10),
    tolerance (1e-10),
    forcing_value (0.5),
    residual_norm (0.0),
    step_norm (0.0),
    iteration (0),
    norm_t (L2),
    t (t_),
    sigma(sigma_),
    theta (theta_),
    thetaPG (thetaPG_),
    theta_min (theta_min_),
    theta_max (theta_max_),
    max_back_it(max_back_it_),
    verbose (verbose_),
    b1(b1_),
    b2(b2_),
    filename ("output.txt")
  {
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
  };

  /// Set nonlinear problem to solve.
  void
  set_problem (abstract_nonlinear_problem *problem_);

  /// Set forcing term.
  void
  set_forcing_term (abstract_forcing_term *forcing_);

  /// Set initial guess.
  void
  set_initial_guess
  (std::vector<double> &initial_guess_);

  /// Set bounds.
  void
  set_bounds
  (std::vector<double> &b1_, std::vector<double> &b2_);

  /// Projection
  void
  projection
 (std::vector<double> &x);

  /// Solve system.
  int
  solve ();

  /// Set maximum number of nonlinear solver iterations.
  void
  set_max_iterations (int max_iter_)
  { max_iter = max_iter_; }

  /// Set tolerance of nonlinear solver.
  void
  set_tolerance (double tolerance_)
  { tolerance = tolerance_; }

  /// Set minimum residual of nonlinear solver.
  void
  set_min_residual (double min_residual_)
  { min_residual = min_residual_; }

  /// Set type of norm used by nonlinear solver.
  void
  set_norm_type (norm_type norm_t_)
  { norm_t = norm_t_; }

  /// Set backtracking parameters.
  /// \f$ t, \sigma \in (0, 1) \f$
  /// \f$ \theta \in [ \theta_{min}, \theta_{max}] \f$
  void
  set_backtracking_parameters
  (double t_,
   double sigma_, 
   double theta_min_,
   double theta_max_)
  {
    t = t_;
    sigma = sigma_;
    theta_min = theta_min_;
    theta_max = theta_max_;
  }

  /// Set theta of the projected Newton direction.
  void
  set_thetaPN
  (double  theta_)
  {
   theta = theta_;
  }
 
  /// Set theta of the projected gradient direction.
  void
  set_thetaPG
  (double  thetaPG_)
  {
    thetaPG = thetaPG_;
  }


  /// Set backtracking max iterations.
  void
  set_backtracking_max_it
  (int  max_back_it_)
  {
    max_back_it = max_back_it_;
  }

  /// \f$ \theta \f$ was chosen to minimize over
  /// \f$ [\theta_{min}, \theta_{max}] \f$ the quadratic
  /// \f$ p (\theta) \f$ for which \f$ p (0) = g (0), p' (0) = g' (0)\f$
  /// and \f$ p (1) = g (1)\f$, where
  /// \f$ g (\theta) = ||F (x_k + \theta * s_k)||_2^2 \f$
  ///
  /// \f$ a = ||F (x_k + s_k)||_2^2 - ||F (x_k)||^2_2 -
  /// 2*F (x_k)^TF' (x_k)s_k \f$
  ///
  /// \f$ b = 2F (x_k)^TF' (x_k)s_k \f$
  ///
  /// \f$ c = ||F (x_k)||^2_2 \f$
  void
  theta_choice (double a, double b, double c);

  /// Set maximum number of linear solver iterations.
  void
  set_max_iterations_of_linear_solver
  (int max_iteration);

  /// Set initial guess of linear solver.
  void
  set_initial_guess_of_linear_solver
  (std::vector<double> &initial_guess);

  /// Set initial tolerance of linear solver.
  void
  set_initial_tolerance_of_linear_solver
  (double initial_tolerance);

  /// Set type of iterative method used by linear solver.
  void
  set_iterative_method_of_linear_solver
  (const std::string &iterative_method);

  /// Set number of iterations for which the iterative method  
  /// used by linear solver is restarted.
  void
  set_restart_iterations_of_linear_solver
  (int restart_iterations);

  /// Set options of iterative method used by linear solver.
  void
  set_options_iterative_method_of_linear_solver
  (const std::string &options_iterative_method);

  /// Set type of preconditioner used by linear solver.
  void
  set_preconditioner_of_linear_solver
  (const std::string &preconditioner);

  /// Set type of convergence condition used by linear solver.
  void
  set_convergence_condition_of_linear_solver
  (const std::string &convergence_condition);

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
  get_result_residual_norm (double &residual_norm_);

  /// Get the solution found by nonlinear solver.
  void
  get_result_solution (std::vector<double> &solution);

  /// Get number of iterations when solve () ends
  void
  get_result_iterations (int &iterations_);

  /// Cleanup memory
  void
  cleanup ();

  /// Set the name of output file.
  void
  set_output_filename (const std::string &filename_);

};

#endif
