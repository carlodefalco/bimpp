/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem:  
   \f[ 
      \begin{cases}
        \partial_t m - \mu \; div  (m \nabla p) = G(p) m \\
        \partial_t n - \nu \; div (n \nabla p) = 0 
      \end{cases}
   \f]
with \f$ p := K_{\gamma}(n+m)^{\gamma} \f$ , \f$ K_{\gamma} := \frac{\gamma + 1}{\gamma} \f$,
\f$ G(p) := \frac{200}{\pi} \arctan( 4 (p - P_M)) \f$ , 
\f$ m \f$ local density of dividing cells (tumor cells), 
\f$ n \f$ local density of non-dividing cells (not tumor cells) 

 \f$ P_M = 30 \f$, \f$ \gamma = 30 \f$.

  Initial conditions: 
  \f$ m(x,y, t=0) := 0.1 e^{-5 \times 10^{-1}(x^2 + y^2)} \f$,
  \f$  n(x,y, t=0) := 0.8 e^{-5 \times 10^{-7}(x^2 + y^2)} \f$.

  Boundary conditions: 
   Neumann homogeneous 

  Exact Solution: 
 
  Linear Solver: lis or mumps

  NonLinear Solver: projected_Newton_method_and_gradient_direction

  Forcing Term: forcing_type3 (1, 2, 0.9)
*/

#include <lis.h>
#include "mumps_class.h"
#include "lis_class.h"
#include "linear_solver.h"
#include "nonlinear_solver.h"
#include "backtracking_inexact_newton_class.h"
#include "adaptive_inexact_newton_class.h"
#include "projected_Newton_method_and_gradient_direction_class.h"
#include "abstract_nonlinear_problem.h"
#include "abstract_forcing_term.h"
#include "forcing_class.h"
#include "bim_config.h"
#include "tumor_growth_class.h"
#include "tmesh.h"
#include "bim_sparse.h"
#include "quad_operators.h"
#include <math.h> 
#include <limits>

#define NUM_CYCLES  6
#define DT  0.01
#define NT  40
#define MIN_RESIDUAL 1e-6
#define MAX_IT 50


constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] = {0., 0., 45., 0.,  45., 45., 0., 45.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = {1, 2, 3, 4, 1};

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_CYCLES; }


int main (int argc, char **argv)
{
  
  MPI_Init (&argc, &argv);
  
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
 
  // linear_solver *lin_solver = new mumps ();

   linear_solver *lin_solver = new lis ();
  
   //  nonlinear_solver *solver =
   // new backtracking_inexact_newton (lin_solver);

   nonlinear_solver *solver =
     new projected_Newton_method_and_gradient_direction (lin_solver);


   
  /// Problem parameters 
  double mu = 1;
  double nu = 2;
  double gamma = 30, PM = 30; 
  double Lx = 0, Rx = 45;
  
  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  int recursive = 1;

  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);
  
  tmesh::idx_t n_nodes = tmsh.num_global_nodes ();
  
  std::vector<double> uold (2 * n_nodes), mold (n_nodes), nold (n_nodes);
  std::vector<double> u (2 *n_nodes), p (n_nodes) ;
  
  /// Implementation of the initial condition of m and n
  double x = 0, y = 0; 
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
	{
	  if (! quadrant->is_hanging (ii))
	    {
	      x = quadrant->p(0, ii);
	      y = quadrant->p(1, ii);
              
	      mold[quadrant->gt(ii)] = 0.1 * std::exp (-0.5* (std::pow(x, 2)
							      + std::pow(y, 2)));             
	      nold[quadrant->gt(ii)] = 0.8 * std::exp (-5.0e-7 * (std::pow(x, 2)
								  + std::pow(y, 2)));             
              
	    }
	}
    }

  for (int i = 0; i < 2 * n_nodes; ++i)
    uold[i] = i < n_nodes ? mold [i] : nold [i-n_nodes];

  for (int i = 0; i <  n_nodes; ++i)
    p[i] = (gamma + 1) / gamma * std::pow (mold[i] + nold[i], gamma);

  
  tmsh.octbin_export ((std::string("tumor_growth_m_")
		       + std::to_string(0)).c_str() , mold );
  tmsh.octbin_export ((std::string("tumor_growth_n_")
			   + std::to_string(0)).c_str() , nold );
  tmsh.octbin_export ((std::string("tumor_growth_p_")
		       + std::to_string(0)).c_str() , p );
  
  /// Time parameters
  double t = 0, T = 1;
  double dt = DT;
  double dt_original = dt;
  int nt = NT;
  std::vector<double> t_save(nt + 1), tstore;

  t_save[0] = t;
  for (int i = 1; i < nt + 1; ++i)
    t_save [i] = t_save[i-1] + (T-t)/nt;
  
  abstract_nonlinear_problem *t_growth = new tumor_growth (mu, nu, t+dt, dt, uold, &tmsh, gamma,
							   PM);
  abstract_forcing_term *forcing = new  forcing_type3 (1, 2, 0.9);

  if (rank == 0)
    {
      std::cout << "\n\n*****\nNon Linear Test\n*****\n";
      
      solver->set_output_filename ("Tumor_growth.txt");
      
      std::cout << "Using solver named "
                << solver->solver_name ()
                << std::endl;
  
    }

  /// Set the parameters for projected Newton with gredient direction
  if (solver->solver_name() == "projected_Newton_method_and_gradient_direction")
    {
      // Bounds of the variables
      std::vector<double> b1, b2;
      b1.assign (2*n_nodes, 0);
      b2.assign (2*n_nodes,  std::numeric_limits<double>::infinity());
      auto tmp = static_cast<projected_Newton_method_and_gradient_direction*> (solver);
      tmp->set_bounds (b1,b2);
          
      // Backtracking parameters
      tmp->set_backtracking_parameters (1e-4, 1e-4, 0.5, 0.8, 0.001, 1);
      tmp->set_backtracking_max_it (15);
    }
      
  /// Set the parameters for backtracking inexact Newton
  if (solver->solver_name () == "Backtracking Inexact Newton")
    static_cast<backtracking_inexact_newton*> (solver)->
      set_backtracking_parameters (1e-4, 0.1, 0.5);
   
  solver->set_problem (t_growth);
  solver->set_forcing_term (forcing);
  solver->set_initial_guess (uold);
    
  
  solver->set_max_iterations (MAX_IT);
  solver->set_tolerance (1e-10);
  solver->set_min_residual (MIN_RESIDUAL);
  solver->set_norm_type (L2);
  
  if (solver->linear_solver_type () == "iterative")
    {
      solver->set_max_iterations_of_linear_solver (2 * n_nodes);
      solver->set_iterative_method_of_linear_solver ("GMRES");
      solver->set_restart_iterations_of_linear_solver (2 * n_nodes);
      solver->set_initial_tolerance_of_linear_solver (0.5);
      solver->set_convergence_condition_of_linear_solver ("norm2_of_rhs");
    }
  
  int half = 0;
  bool converged;
  double residual_norm;
  int nonlinear_it; 
  double dtold = dt;


  t += dt;

  for (int its = 1; its < nt + 1; ++its) // nt + 1
    {
      while (t < t_save[its] )
	{
	  if (rank == 0)
	    std::cout << "############## TIME ############## : "<< t << std::endl;
	  
	  dt = dt_original;
	  if (t + dt > t_save[its] )
	    dt = t_save[its] - t;
	    

	  //// NON MI PIACE /////////////////////////////////
	  
	  /// Set the parameters for projected Newton with gredient direction
	  if (solver->solver_name() == "projected_Newton_method_and_gradient_direction")
	    static_cast<projected_Newton_method_and_gradient_direction*> (solver)
	      ->set_backtracking_parameters (1e-4, 1e-4, 0.5, 0.8, 0.001, 1);
	      
	  /// Set the parameters for backtracking inexact Newton
	  if (solver->solver_name () == "Backtracking Inexact Newton")
	    static_cast<backtracking_inexact_newton*> (solver)->
	      set_backtracking_parameters (1e-4, 0.1, 0.5);
	    
	  ////////////////////////////////////////////////

	  static_cast<tumor_growth*>(t_growth)->set_matrices_structure (); 

	  converged = solver->solve ();
	  
	  solver->get_result_residual_norm (residual_norm);
              
	  if (rank == 0 && converged == 0)
	    {
	      std::cout << "The algorithm did't converge"<< std::endl;
	      exit(-1);
	    }
	  solver->get_result_solution (u);
	  double min_u = 10;
	  for (int i = 0; i < u.size(); ++i)
	    min_u = std::min (u[i], min_u);
	  /* if (min_u < 0 || residual_norm > MIN_RESIDUAL)
	    {
	      half = half + 1; 
	      dt = dtold / 2;
	      t = t - dtold + dt ;
	      static_cast<tumor_growth*>(t_growth)->set_t_dt (t, dt);
	    }
	  else
	  {*/
	      for (int i = 0; i < 2 * n_nodes; ++i)
		uold[i] = u[i];
	      solver->set_initial_guess (uold);
	      tstore.push_back (t);
	      t += dt;
	      static_cast<tumor_growth*>(t_growth)->set_initial_condition(uold);
	      static_cast<tumor_growth*>(t_growth)->set_t_dt (t, dt);
	      //  if (rank == 0)
	      //	std::cout << "Times in which dt was halved : " << half <<std::endl;
	      // half = 0;
	      // }
	  dtold = dt;  

	}
      
      for (int i = 0; i < n_nodes; ++i)
	{
	  mold[i] = uold[i];
	  nold[i] = uold[i + n_nodes];
	  p[i] = (gamma + 1) / gamma * std::pow (mold[i] + nold[i], gamma);
	}
      tmsh.octbin_export ((std::string("tumor_growth_m_")
			   + std::to_string(its)).c_str() , mold );
      tmsh.octbin_export ((std::string("tumor_growth_n_")
			   + std::to_string(its)).c_str() , nold );
      tmsh.octbin_export ((std::string("tumor_growth_p_")
			   + std::to_string(its)).c_str() , p);
         
    }

  solver->cleanup ();
  
  
  MPI_Finalize ();
}


