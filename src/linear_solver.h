/*
  Copyright (C) 2015 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file linear_solver.h
  \brief Generic interface for a linear solver.
*/

#ifndef HAVE_LINEAR_SOLVER_H
#define HAVE_LINEAR_SOLVER_H 1

#include <string>

/// Generic interface for a linear solver.
class linear_solver
{
private :

  const std::string name, type;

protected :

  linear_solver (const char *name_, const char *type_) :
    name (name_), type (type_) { };

public :

  /// Format for sparse matrix structure.
  enum matrix_format_t {aij = 0, csr = 1};

  /// \brief Set-up the matrix structure.
  /// \details Memory for lhs and rhs is allocated on the
  /// master (rank == 0) ans slave (rank != 0) nodes.
  /// The nonzero pattern for the lhs is set, but matrix
  /// entries are not yet assigned.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  set_lhs_structure
  (int number_of_rows,
   std::vector<int> &i_rows,
   std::vector<int> &j_columns,
   matrix_format_t format = aij) = 0;

  /// \brief Perform analysis steps required prior to factorization
  /// (e.g. reordering, partitioning, etc.).
  /// \details Return values are implementation dependent,
  /// but a positive return value always indicates success.
  /// The lhs structure must be set before this step, but
  /// lhs entries are not necessarily required.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual int
  analyze () { return 1; };

  /// \brief Set values of the lhs matrix entries.
  /// \details Must be called on the master (rank == 0)
  /// node only.
  virtual void
  set_lhs_data
  (std::vector<double> &matrix_entries) = 0;

  /// \brief Perform matrix factorization if required.
  /// \details Return values are implementation dependent,
  /// but a positive return value always indicates success.
  /// The lhs structure and values must be set before this step.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual int
  factorize () { return 1; };

  /// \brief Set the rhs vector value.
  /// \details rhs must not be freed after set_rhs
  /// as it will be used internally to store
  /// the system solution. The whole vector
  /// must be allocated on the master (rank == 0)
  /// node.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  set_rhs (std::vector<double> &rhs) = 0;

  /// \brief Solve the system
  /// \details rhs must not be freed after solve
  /// as it will be used internally to store
  /// the system solution.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual int
  solve () = 0;

  /// \brief Cleanup memory.
  /// \details Must be called on the master (rank == 0)
  /// and slave (rank != 0) nodes at the same time.
  /// After invoking this method the object should not be used anymore.
  virtual void
  cleanup () { };

  /// \brief Set the initial guess of linear solver.
  /// \details Must be called only if linaer solver is an iterative method.
  /// Must be called on the master (rank == 0) node only.
  virtual void
  set_initial_guess (std::vector<double> &initial_guess) { };

  /// \brief Set maximum number of iterations of linear solver.
  /// \details Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_max_iterations (int max_iter) { };

  /// \brief Get maximum number of iterations of linear solver.
  /// \details Must be called only if linear solver is an iterative method.
  virtual void
  get_max_iterations (int &max_iter) { };

  /// \brief Set tolerance of linear solver.
  /// \details Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual void
  set_tolerance (double tolerance) { };

  /// \brief Get tolerance of linear solver.
  /// \details Must be called only if linear solver is an iterative method.
  virtual void
  get_tolerance (double &tolerance) { };

  /// \brief Set type of iterative method used to find solution of system.
  /// \details Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  ///
  /// Possible input are:
  ///
  /// "Conjugate Gradient" for Conjugate Gradient Method
  ///
  /// "Biconjugate Gradient" for BiConjugate Gradient Method
  ///
  /// "Bicg Stabilized" for BiConjugate Gradient Stabilized Method
  ///
  /// "Jacobi" for Jacobi Method
  ///
  /// "Gauss Seidel" for Gauss Seidel Method
  ///
  /// "SOR" for SOR Method
  ///
  /// If iterative method sent to solver is invalid,
  /// the program will run with default
  /// iterative method "Biconjugate Gradient".
  virtual void
  set_iterative_method
  (const std::string &type_of_iterative_method) { };

  /// \brief Get type of iterative method sent to solver.
  /// \details Must be called only if linear solver is an iterative method.
  virtual void
  get_iterative_method
  (std::string &type_of_iterative_method) { };

  /// \brief Set preconditioner of linear solver.
  /// \details Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  ///
  /// Possible input are:
  ///
  /// "none" for solve system without preconditioner.
  ///
  /// "jacobi" for Jacobi preconditioner.
  ///
  /// "ssor" for Symmetric Successive Over-Relaxation.
  ///
  /// If preconditioner sent to solver is invalid,
  /// the program will run without preconditioner.
  virtual void
  set_preconditioner (const std::string &preconditioner) { };

  /// \brief Get preconditioner sent to solver.
  /// \details Must be called only if linear solver is an iterative method.
  virtual void
  get_preconditioner (std::string &preconditioner) { };

  /// \brief Set convergence condition of iterative method.
  /// \details Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  /// Possible input are:
  ///
  /// "norm2_of_residual" for \f$ ||b-Ax||_2 <= tol * ||b-Ax_0||_2 \f$
  ///
  /// "norm2_of_rhs" for \f$ ||b-Ax||_2 <= tol * ||b||_2 \f$
  ///
  /// If convergence condition sent to solver is invalid,
  /// the program will run with default convergence condition
  /// "norm2_of_residual".
  virtual void
  set_convergence_condition
  (const std::string &convergence_condition) { };

  /// \brief Get convergence condition of iterative method sent to solver.
  /// \details Must be called only if linear solver is an iterative method.
  virtual void
  get_convergence_condition
  (std::string &convergence_condition) { };

  /// Return the name of the specific implementation.
  const std::string&
  solver_name () { return name; }

  /// \brief Return the type (either "iterative" or "direct")
  /// of the specific implementation.
  const std::string&
  solver_type () { return type; }

  /// \brief Return the preferred base used for indexing
  /// \details usually libraries written in C such as LIS
  /// prefer 0-based indexing while libraries written
  /// in fortran prefer 1-based indexing.
  virtual int
  get_index_base ()
  { return 0; }

};

#endif
