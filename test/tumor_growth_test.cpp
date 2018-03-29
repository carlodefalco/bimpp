/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem:  
  \f[ -div (|\nabla (u)|^{p-2} \nabla (u)) = f \f]

  \f[ u = 1/q \cdot (0.5^q -
     [ (x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2)]^{q/2} \:on \:boundary \f]

  \f[ f = 3.0 \f]

  \f[ p = 3.0 \f]

  Exact Solution: 
  \f[  u = 1/q \cdot (0.5^q -
           [ (x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2)]^{q/2} \f]

  Linear Solver: lis

  NonLinear Solver: backtracking_inexact_newton

  Forcing Term: forcing_type3 (1, 2, 0.9)
*/

#include <lis.h>
#include "mumps_class.h"
#include "lis_class.h"
#include "linear_solver.h"
#include "operators.h"
#include "mesh.h"
#include "nonlinear_solver.h"
#include "backtracking_inexact_newton_class.h"
#include "projected_Newton_method_and_gradient_direction_class.h"
#include "backtracking_inexact_newton_example8_class.h"
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

constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] = {0., 0., 45., 0.,  45., 45., 0., 45.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = {1, 2, 3, 4, 1};

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return 1; }

void
run_test_problem (nonlinear_solver *solver);

int main (int argc, char **argv)
{
  
  MPI_Init (&argc, &argv);
  
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
 
  /// Limits
  std::vector<double> b1, b2;
  b1.assign (50 , 0);
  b2.assign (50 ,  std::numeric_limits<double>::infinity());

  linear_solver *lin_solver = new lis ();
  //  linear_solver *lin_solver = new mumps ();
   
  nonlinear_solver *solver =
    new projected_Newton_method_and_gradient_direction (lin_solver, b1 , b2 );
  
  run_test_problem (solver);
  
  MPI_Finalize ();
}

void
run_test_problem (nonlinear_solver *solver)
{
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  
  //std::vector<double> exactsolution;

  //  std::vector<int>    bnodes;
  //std::vector<double> vnodes;
  
  // std::vector<double> uold;

  //std::vector<double> fcoeff;

  /// Problem's parameters 
  double mu = 1;
  double nu = 2;
  double Lx = 0, Rx = 45;  
  int number_cycles = 2; 

  /// Generate the mesh in 2D
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  tmsh.set_replace_fun (tmesh::userint_replace);
  int recursive , partforcoarsen, balence;
  recursive = 0;
  partforcoarsen = 1;
  
  for (int cycle = 0; cycle < number_cycles; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  tmsh.vtk_export ("mesh_tumor_growth");
  tmsh.save("file_mesh_tumor_growth");
  int n_nodes = tmsh.num_global_nodes ();
  
  std::vector<double> tstore, uold (2 * n_nodes), mold (n_nodes), nold (n_nodes), u (2*n_nodes) ;

  
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
              
	      mold[quadrant->gt(ii)] = 0.1 * std::exp (-5 * std::pow (10, -1) 
						       * (std::pow(x, 2) + std::pow(y, 2)));             
	      nold[quadrant->gt(ii)] = 0.8 * std::exp (-5 * std::pow (10, -7) 
						       * (std::pow(x, 2) + std::pow(y, 2)));             
              
	    }
	}
    }
  if (rank == 0)
    {
      for (int i = 0; i < 2 * n_nodes; ++i)
	{
	  uold[i] = i < n_nodes ? mold [i] : nold [i-n_nodes];
	}
    }

  if (rank == 0)
    {
      tmsh.octbin_export ((std::string("tumor_growth_m_")
			   + std::to_string(0)).c_str() , mold );
      tmsh.octbin_export ((std::string("tumor_growth_n_")
			   + std::to_string(0)).c_str() , nold );
    }
  
  /// Time's parameters
  double t = 0, T = 1;
  double dt = 0.1;
  double dt_original = dt;
  int nt = 2;
  std::vector<double> t_save(nt + 1);
  if (rank == 0)
    {
      t_save[0] = t;
      for (int i = 1; i < nt + 1; ++i)
	{
	  t_save [i] = t_save[i-1] + (T-t)/nt;
	}
    }
  abstract_nonlinear_problem *t_growth = new tumor_growth (mu, nu, t+dt, dt, uold);
  abstract_forcing_term *forcing = new  forcing_type3 (1, 2, 0.9);

  if (rank == 0)
    {
      std::cout << "\n\n*****\nNon Linear Test\n*****\n";
      
      solver->set_output_filename ("Tumor_growth.txt");

      std::cout << "Using solver named "
                << solver->solver_name ()
                << std::endl;

      t_growth->read_mesh ("file_mesh_tumor_growth");
    }

  if (rank == 0)
    {
      if (solver->solver_name () == "Backtracking Inexact Newton")
        ((backtracking_inexact_newton *) solver)->
          set_backtracking_parameters (1e-4, 0.1, 0.5);
      solver->set_problem (t_growth);
      solver->set_forcing_term (forcing);
      solver->set_initial_guess (uold);
    }

  solver->set_max_iterations (20);
  solver->set_tolerance (1e-10);
  solver->set_min_residual (1e-10);
  solver->set_norm_type (L2);
  
  if (solver->linear_solver_type () == "iterative")
    {
      solver->set_max_iterations_of_linear_solver (1000);
      solver->set_iterative_method_of_linear_solver
	("Conjugate Gradient");
      solver->set_initial_tolerance_of_linear_solver (0.5);
      solver->set_convergence_condition_of_linear_solver ("norm2_of_rhs");
    }

  ///////////////////////////////////////////////////////////////////////////////////////////
  
  bool converged;
  t += dt;
  for (int its = 1; its < nt + 1; ++its)
    {
      while (t < t_save[its] )
	{
	  if (rank == 0)
	    {
	      std::cout << "############## TIME ############## : "<< t << std::endl;
	      dt = dt_original;

	      if (t + dt > t_save[its] )
		dt = t_save[its] - t;
	    }
	  
	  converged = solver->solve ();
	  if ( rank == 0)
	    {
	      if (converged == 0)
		{
		  std::cout << "The algorithm did't converge"<< std::endl;
		  exit(-1);
		}
	      solver->get_result_solution (u);
	  
	      for (int i = 0; i < 2 * n_nodes; ++i)
		uold[i] = u[i];
	      
	      solver->set_initial_guess (uold);

	      t += dt;
	      ((tumor_growth*) t_growth)->set_initial_condition(uold);
	      ((tumor_growth*) t_growth)->set_t_dt (t, dt);
	    }
	}
      if ( rank == 0)
	{
	  for (int i = 0; i < n_nodes; ++i)
	    {
	      mold[i] = uold[i];
	      nold[i] = uold[i + n_nodes];
	    }
	  tmsh.octbin_export ((std::string("tumor_growth_m_")
			       + std::to_string(its)).c_str() , mold );
	  tmsh.octbin_export ((std::string("tumor_growth_n_")
			       + std::to_string(its)).c_str() , nold );
	}
    }

  ///////////////////////////////////////////////////////////////////////////////////////////

  solver->cleanup ();
  
}
