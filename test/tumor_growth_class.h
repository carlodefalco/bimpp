/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file plaplacian_class.h
  \brief Interface for nonlinear problem
  \f$ -div (|\nabla u|^{p-2}\nabla u) = f \f$
*/

#ifndef HAVE_PLAPLACIAN_CLASS_H
#define HAVE_PLAPLACIAN_CLASS_H

#include "abstract_nonlinear_problem.h"
#include "operators.h"
#include "bim_sparse.h"
#include "tmesh.h"
#include <math.h>
#include "quad_operators.h"


/// \brief Interface for nonlinear problem
/// \f$ div (|\nabla u|^{p-2}\nabla u) = f \f$
class tumor_growth : public abstract_nonlinear_problem
{
private :

  /// Mobilities.
  double mu;
  double nu;
  
  /// gamma.
  double g; 
  
  /// Homeostatic pressure.
  double PM;
  
  /// Current time 
  double t;

  /// Time step
  double dt;

  /// Previous solution
  std::vector<double> uold;
  
  /// Stores the exact solution of the nonlinear problem.
  std::vector<double> exact_solution;

  /// Stores rhs values of nonlinear problem.
  std::vector<double> f;

  /// Stores boundary values.
  std::vector<double> boundary_values;

  /// Stores boundary nodes.
  std::vector<int> boundary_nodes;

  /// mesh 
  tmesh tmsh;  // VA GESTITO COSÌ ? 

public :

  /// Default costructor.
  tumor_growth (double mu_, double nu_, double t_, double dt_, std::vector<double> &uold_,
		double g_=30, double PM_ = 30) :
    abstract_nonlinear_problem ("tumor_growth"),
      mu (mu_), nu(nu_), t(t_), dt(dt_), uold(uold_), g(g_), PM(PM_) { };

  /// Read mesh.
  void
  read_mesh (const std::string &mesh_name);

  /// Set t and dt
  void
  set_t_dt (const double &t_, const double &dt_)
  {
    t = t_ ;
    dt = dt_;
  }
  
  /// Set the exact solution of nonlinear problem.
  void
  set_exact_solution (const std::vector<double> & exact_solution_);

  /// Set the solution at the previous time step.
  void
  set_initial_condition (const std::vector<double> & uold_);

  /// Set rhs values of nonlinear problem.
  void
  set_rhs_values (const std::vector<double> & f_);

  /// Set boundary values and nodes.
  void
  set_boundary_conditions
    (std::vector<double> &boundary_values_,
     std::vector<int> &boundary_nodes_);

  /// Compute lhs and rhs of linearized nonlinear problem in guess.
  void
  operator () (sparse_matrix& lhs,
               std::vector<double>& rhs,
               const std::vector<double>& guess);
 
  /// Valued nonlinear functional in guess.
  void
  operator () (std::vector<double>& functional,
               const std::vector<double>& guess);

  /// Get the exact solution of nonlinaer problem.
  void
  get_exact_solution (std::vector<double> &exact_solution_);
};

#endif
