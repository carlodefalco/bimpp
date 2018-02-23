/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem:
  \f[ x_1^2 - 1 = 0 \f]
  \f[ x_1 - x_2^3 = 0\f]
  \f[ ... \f]
  \f[ x_{n-2} - x_{n-1}^3 = 0 \f]
  \f[ x_{n-1} - x_{n} = 0 \f]
  \f[ x \in \Omega \f]
  
  Exact Solution:
  \f[ x = (1, ... , 1)^T \f]
  
*/

#ifndef HAVE_PLAPLACIAN_CLASS_H
#define HAVE_PLAPLACIAN_CLASS_H

#include "abstract_nonlinear_problem.h"
#include "operators.h"
#include "bim_sparse.h"

///  \brief Interface for nonlinear problem is 
///   taken from Example 8 in the paper
///  "Globalization technique for ptojeced Newton-Krylov methods"

class example8 : public abstract_nonlinear_problem
{
private :


  /// Stores the exact solution of the nonlinear problem.
  std::vector<double> exact_solution;

  /// Stores rhs values of nonlinear problem.
  std::vector<double> f;

public :
  /// Dimention of the nonlinear system.
  unsigned int n ; 

  /// Default costructor.
  example8 (unsigned int n_) :
    abstract_nonlinear_problem ("example8"),
    n (n_) { };

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
  operator () (std::vector<double>& functional,
               const std::vector<double>& guess);

  /// Get the exact solution of nonlinear problem.
  void
  get_exact_solution (std::vector<double> &exact_solution_);
  

};

#endif
