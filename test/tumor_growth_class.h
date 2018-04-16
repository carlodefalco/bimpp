/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file tumor_growth_class.h
  
  \brief Interface for nonlinear problem
   \f[ 
      \begin{cases}
        \partial_t m - \mu \; div  (m \nabla p) = G(p) m \\
        \partial_t n - \nu \; div (n \nabla p) = 0 
      \end{cases}
   \f]
with \f$ p := K_{\gamma}(n+m)^{\gamma} \f$ , \f$ K_{\gamma} := \frac{\gamma + 1}{\gamma} \f$,

\f$ G(p) := \frac{200}{\pi} \arctan( 4 (p - P_M)) \f$ , 

\f$ m \f$ local density of dividing cells (tumor cells), 
\f$ n \f$ local density of non-dividing cells (not tumor cells) .

*/

#ifndef HAVE_TUMOR_GROWTH_H
#define HAVE_TUMOR_GROWTH_H 1

#include "abstract_nonlinear_problem.h"
#include "bim_sparse.h"
#include "tmesh.h"
#include <math.h>
#include "quad_operators.h"


/// \brief Interface for nonlinear problem
///   is about the growth of tumor cells
/// \details It's taken form the paper "On interfaces between cell populations with different 
///  mobilities" of Lorenzi, Lorz and Perthame.
/// \f[ 
///      \begin{cases}
///        \partial_t m - \mu \; div  (m \nabla p) = G(p) m \\
///        \partial_t n - \nu \; div (n \nabla p) = 0 
///      \end{cases}
///   \f]
/// with \f$ p := K_{\gamma}(n+m)^{\gamma} \f$ , \f$ K_{\gamma} := \frac{\gamma + 1}{\gamma} \f$,
///
/// \f$ G(p) := \frac{200}{\pi} \arctan( 4 (p - P_M)) \f$ ,
/// 
/// \f$ m \f$ local density of dividing cells (tumor cells), 
/// \f$ n \f$ local density of non-dividing cells (not tumor cells) .

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
  tmesh *tmsh;
  
  tmesh::idx_t n_nodes ; 
  tmesh::idx_t n_elements ;

  /// Matrixes and vectors used 
  sparse_matrix Amm, Ann, mass, Sm, Sn, mat_temp, M;
  std::vector<double> m, n, mold, nold;
  std::vector<double> diffm, diffn, p;
  std::vector<double> G, mdGdm, mdGdn;
  std::vector<double> tempm , tempn, tempb;
  std::vector<double> ecoeff;
  std::vector<double> ncoeff;

public :

  /// Default costructor.
  tumor_growth (double mu_, double nu_, double t_, double dt_, std::vector<double> &uold_,
	        tmesh *tmsh_, double g_=30, double PM_ = 30 ) :
  abstract_nonlinear_problem ("tumor_growth"),
    mu (mu_), nu (nu_), t (t_), dt (dt_), uold (uold_), tmsh (tmsh_), g (g_), PM (PM_) { };
  
  /// Read mesh.
  void
    read_mesh (const std::string &mesh_name) {};

  /// Set t and dt
  void
  set_t_dt (const double &t_, const double &dt_)
  {
    t = t_ ;
    dt = dt_;
  }

  /// Initialize the structure of the sparse matrices used in operator()
  void 
  set_matrices_structure () ;

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

  /// Set exact solution of the nonlinear problem.
  /// Must be called on the master (rank == 0) node only.
  void
  set_exact_solution (const std::vector<double> exact_solution) { };

  /// Get the exact solution of the nonlinear problem.
  /// Must be called on the master (rank == 0) node only.
  void
  get_exact_solution (std::vector<double> &exact_solution) { };

  
};

#endif
