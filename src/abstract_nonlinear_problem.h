/*
  Copyright (C) 2015 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file abstract_nonlinear_problem.h
  \brief Generic interface for a nonlinear problem
*/

#ifndef HAVE_ABSTRACT_NONLINEAR_PROBLEM_H
#define HAVE_ABSTRACT_NONLINEAR_PROBLEM_H

#include <string>
#include "mesh.h"
#include "bim_sparse.h"

/// Generic interface for a nonlinear problem
class abstract_nonlinear_problem
{
private :

  /// The name of specific nonlinear problem.
  const std::string name;

protected :

  abstract_nonlinear_problem (const char *name_) :
    name (name_) { };

public :

  /// Read mesh used for the nonlinear problem.
  /// Must be called on the master (rank == 0) node only.
  virtual  void
  read_mesh (const std::string &mesh_name) { };

  /// Set exact solution of the nonlinear problem.
  /// Must be called on the master (rank == 0) node only.
  virtual void
  set_exact_solution (const std::vector<double> exact_solution) { };

  /// Set values of the function on rhs
  /// of the nonlinear problem \f$ A(x) = f \f$.
  /// Must be called on the master (rank == 0) node only.
  virtual void
  set_rhs_values (const std::vector<double> & f_) { };

  /// Set values and nodes of boundary conditions
  /// of the nonlinear problem.
  /// Must be called on the master (rank == 0) node only.
  virtual void
  set_boundary_conditions
  (std::vector<double> &boundary_values,
   std::vector<int> &boundary_nodes) { };

  /// Operator that computes lhs and rhs
  /// of the linearized problem \f$ DF(x) = - F(x)\f$ valued in guess.
  /// Must be called on the master (rank == 0) node only.
  virtual void
  operator () (sparse_matrix &lhs,
               std::vector<double> &rhs,
               const std::vector<double> &guess) = 0;


  /// Operator that computes \f$ F(x) = A(x) - f \f$ valued in guess.
  /// Must be called on the master (rank == 0) node only.
  virtual void
  operator () (std::vector<double> &functional,
               const std::vector<double> &guess) = 0;

  /// Get the exact solution of the nonlinear problem.
  /// Must be called on the master (rank == 0) node only.
  virtual void
  get_exact_solution (std::vector<double> &exact_solution) { };

  /// Returns the name of the nonlinear problem.
  const std::string&
  problem_name () { return name; }

  /// Mesh used in the nonlinear problem.
  mesh msh;
};

#endif
