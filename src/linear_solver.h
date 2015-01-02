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

class linear_solver
{
public :

  /// Set-up the matrix structure.
  /// Memory for lhs and rhs is allocated on the
  /// master (rank == 0) ans slave (rank <> 0) nodes.
  /// The nonzero pattern for the lhs is set, but matrix
  /// entries are not yet assigned.
  /// Must be called on the master (rank == 0)
  /// node only.
  virtual void 
  set_lhs_structure
  (int number_of_rows, 
   std::vector<int> &i_rows, 
   std::vector<int> &j_columns) = 0;

  /// Perform analysis steps required prior to factorization
  /// (e.g. reordering, partitioning, etc.).
  /// Return values are implementation dependent,
  /// but a positive return value always indicates success.
  /// The lhs structure must be set before this step, but
  /// lhs entries are not necessarily required.
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
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
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
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
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
  /// nodes at the same time.
  virtual int
  solve () = 0;

  /// Cleanup memory.
  /// Must be called on the master (rank == 0) and slave (rank <> 0)
  /// nodes at the same time.
  /// After invoking this method the object should not be used anymore.
  virtual void 
  cleanup () { };

};

#endif
