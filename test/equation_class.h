/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file equation_class.h
  \brief interface for nonlinear problem \f$ x^2 = p \f$
 */
#ifndef HAVE_EQUATION_CLASS_H
#define HAVE_EQUATION_CLASS_H

#include "abstract_nonlinear_problem.h"

class equation : public abstract_nonlinear_problem
{
private:

  double value;
  std::vector<double> exact_solution;

public :

  equation (double value_) :
    abstract_nonlinear_problem ("2nd_order_equation"),
    value (value_) { };

  void
  operator () (sparse_matrix& lhs,
               std::vector<double>& rhs,
               const std::vector<double>& guess);
  void
  operator () (std::vector<double>& functiornal,
               const std::vector<double>& guess);
  void
  get_solution (std::vector<double>& sol);

};

#endif

