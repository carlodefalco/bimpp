/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file nonlinear_solver.h
  \brief generic interface for a nonlinear solver.
*/

#ifndef HAVE_NL_SOLVER_H
#define HAVE_NL_SOLVER_H 1

#include <string>
#include "abstract_nonlinear_problem.h"
#include "abstract_forcing_term.h"

class nonlinear_solver
{
private :

  const std::string name;

protected :

  nonlinear_solver (const char *name_) :
    name (name_) { };

public :

  /// Set nonlinear problem to solve.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  set_problem (abstract_nonlinear_problem *problem_) = 0;

  /// Set forcing term which calculates forcing value.
  /// after every iteration of nonlinear solver.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  set_forcing_term (abstract_forcing_term *forcing_term) = 0;

  /// Set the initial guess
  /// to compute first linearization of nonlinear problem.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  set_initial_guess
  (std::vector<double> &initial_guess) = 0;

  /// Solve nonlinear problem.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  /// Return 0 if maximum number of iterations reached
  /// Return 1 if the method converged
  virtual int
  solve () = 0;

  /// Set maximum number of nonlinear solver iterations.
  /// Must be called on master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_max_iterations (int max_iter) = 0;

  /// Set tolerance of nonlinear solver.
  /// nonlinear iterations stop if
  /// \f$ \|x_{new} - x_{old}\| < tolerance \f$
  /// Must be called on master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_tolerance (double tolerance) = 0;

  /// Set minimum residual of nonlinear solver.
  /// Nonlinear iterations stop if
  /// \f$ \|F(x)\| < min_residual \f$
  /// where F is the nonlinear function.
  /// Must be called on master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_min_residual (double min_res) = 0;

  /// Set type of norm used by nonlinear solver.
  /// Must be called on master (rank == 0)
  /// node only.
  virtual void
  set_norm_type (norm_type type_of_norm) = 0;

  /// Set maximum number of linear solver iterations.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_max_iterations_of_linear_solver
  (int max_iteration) { };

  /// Set initial guess of iterative method used by linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) node only.
  virtual void
  set_initial_guess_of_linear_solver
  (std::vector<double> &initial_guess) { };

  /// Set initial tolerance of iterative method used by linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_initial_tolerance_of_linear_solver
  (double initial_tolerance) { };

  /// Set type of iterative method used by linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_iterative_method_of_linear_solver
  (const std::string &iterative_method) { };

  /// Set type of preconditioner used by linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_preconditioner_of_linear_solver
  (const std::string &preconditioner) { };

  /// Set type of convergence condition used by linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_convergence_condition_of_linear_solver
  (const std::string &convergence_condition) { };

  /// Return the name of linear solver used.
  virtual const std::string&
  linear_solver_name () { };

  /// Return the type of linear solver used.
  virtual const std::string&
  linear_solver_type () { };

  /// Get the norm of the solution's residual.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  get_result_residual_norm (double &residual_norm_) { };

  /// Get the solution calculated with nonlinear solver.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  get_result_solution (std::vector<double> &solution) { };

  /// Get number of iterations when solve () ends
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  get_result_iterations (int &iterations_) { };

  /// Cleanup memory.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  /// After invoking this method the object should not be used anymore.
  virtual void
  cleanup () { };

  /// Return the name of the specific implementation.
  const std::string&
  solver_name () { return name; }

  /// Set the name of output file.
  virtual void
  set_output_filename (const std::string &filename) { };

};

#endif
