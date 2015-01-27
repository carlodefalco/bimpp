/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file equation_class.h
  \brief Interface for nonlinear problem \f$ x^2 = p \f$
 */
#ifndef HAVE_EQUATION_CLASS_H
#define HAVE_EQUATION_CLASS_H

#include "abstract_nonlinear_problem.h"

/// Interface for nonlinear problem \f$ x^2 = p \f$
class equation : public abstract_nonlinear_problem
{
private:

  /// rhs of nonlinear equation.
  double value;

  /// Stores the exact solution of equation.
  std::vector<double> exact_solution;

public :

  /// Default costructor.
  equation (double value_) :
    abstract_nonlinear_problem ("2nd_order_equation"),
    value (value_) { };

  /// Set the exact solution of nonlinear problem.
  void
  set_exact_solution (const std::vector<double> exact_solution_);

  /// Compute lhs and rhs of linearized nonlinear problem in guess.
  void
  operator () (sparse_matrix& lhs,
               std::vector<double>& rhs,
               const std::vector<double>& guess);

  /// Valued nonlinear functional in guess.
  void
  operator () (std::vector<double>& functiornal,
               const std::vector<double>& guess);

  /// Get the exact solution of nonlinaer problem.
  void
  get_exact_solution (std::vector<double>& sol);

};

#endif

