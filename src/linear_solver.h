/*
  Copyright (C) 2015 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/
/*! \file linear_solver.h
  \brief generic interface for a linear solver.
*/

#ifndef HAVE_LINEAR_SOLVER_H
#define HAVE_LINEAR_SOLVER_H 1

#include <string>

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

  /// Set-up the matrix structure.
  /// Memory for lhs and rhs is allocated on the
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

  /// Perform analysis steps required prior to factorization
  /// (e.g. reordering, partitioning, etc.).
  /// Return values are implementation dependent,
  /// but a positive return value always indicates success.
  /// The lhs structure must be set before this step, but
  /// lhs entries are not necessarily required.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual int
  analyze () { return 1; };

  /// Set values of the lhs matrix entries.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  set_lhs_data
  (std::vector<double> &matrix_entries) = 0;

  /// Perform matrix factorization if required.
  /// Return values are implementation dependent,
  /// but a positive return value always indicates success.
  /// The lhs structure and values must be set before this step.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual int
  factorize () { return 1; };

  /// Set the rhs vector value.
  /// rhs must not be freed after set_rhs
  /// as it will be used internally to store
  /// the system solution. The whole vector
  /// must be allocated on the master (rank == 0)
  /// node.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void
  set_rhs (std::vector<double> &rhs) = 0;

  /// Set the rhs vector value.
  /// rhs must not be freed after set_rhs
  // as it will be used internally to store
  /// the system solution.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  virtual int
  solve () = 0;

  /// Cleanup memory.
  /// Must be called on the master (rank == 0) and slave (rank != 0)
  /// nodes at the same time.
  /// After invoking this method the object should not be used anymore.
  virtual void
  cleanup () { };

  /// Set max iteration of linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
  /// nodes at the same time.
  virtual void
  set_max_iterations (int max_iter) { };

  /// Get max iteration of linear solver.
  /// Must be called only if linear solver is an iterative method.
  virtual void
  get_max_iterations (int &max_iter) { };

  /// Set tolerance of linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
  /// nodes at the same time.
  virtual void
  set_tolerance (double tolerance) { };

  /// Get tolerance of linear solver.
  /// Must be called only if linear solver is an iterative method.
  virtual void
  get_tolerance (double &tolerance) { };

  /// Set type of linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
  /// nodes at the same time.
  virtual void
  set_linear_solver (const std::string &linear_solver_type) { };

  /// Get type of linear solver.
  /// Must be called only if linear solver is an iterative method.
  virtual void
  get_linear_solver (std::string &linear_solver_type) { };
  

  /// Set preconditioner of linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
  /// nodes at the same time.
  virtual void
  set_preconditioner (const std::string &preconditioner) { };

  /// Get preconditioner of linear solver.
  /// Must be called only if linear solver is an iterative method.
  virtual void
  get_preconditioner (std::string &preconditioner) { };

  /// Set other options of linear solver.
  /// Must be called only if linear solver is an iterative method.
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
  /// nodes at the same time.
  virtual void
  set_other_options (const std::string &other_options) { };

  /// Get other options of linear solver.
  /// Must be called only if linear solver is an iterative method.
  virtual void
  get_other_options (std::string &other_options) { };

  /// Return the name of the specific implementation.
  const std::string&
  solver_name () { return name; }

  /// Return the type (either "iterative" or "direct")
  /// of the specific implementation.
  const std::string&
  solver_type () { return type; }

};

#endif
