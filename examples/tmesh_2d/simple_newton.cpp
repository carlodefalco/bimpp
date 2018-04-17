/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <cmath>

#include <bim_timing.h>
#include <mumps_class.h>
#include <tmesh.h>
#include <quad_operators.h>


constexpr int NUM_REFINEMENTS = 5;
constexpr double MIN_RESIDUAL = 1.e-6;
constexpr int NT              = 10;
constexpr int MAX_IT          = 50;
constexpr double DT           = 1.0e-2;

constexpr double T0    = 0.;
constexpr double T     = 10.;
constexpr double eps_u = 1.0e-2;

constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {0., 0., 1., 0.,  1., 1., 0., 1.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {1, 2, 3, 4, 1};

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

bool
all_non_negative (std::vector<double> vect)
{
  for (int i = 0; i < vect.size () ; ++i)
    if (vect[i] < 0)
      return 0;
  return 1; 
}


int
main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);

  tmesh::idx_t n_nodes = tmsh.num_global_nodes ();
  tmesh::idx_t n_elements = tmsh.num_local_quadrants ();

  /// Time stepping variables
  double t           = T0;
  double told        = T0;
  double tvold       = T0;
  double dt          = DT;
  int nt             = NT;

  std::array<double, NT+1> t_save {t};
  std::vector<double> t_vect;

  auto its = [&t] ()
    {double out = t; t += ((T - T0) / NT); return out; };
  std::generate (t_save.begin (), t_save.end (), its);
  t = T0;
  
  /// Newton variables
  double              residual_norm;
  int                 it_nonlin;
  auto compute_norm = [&residual_norm] (double x)
    { residual_norm += std::pow (x, 2); };

  /// PDE arrays
  // (1/dt) * mass * (u - uold) + laplacian (eps) * u + mass * (exp (u)) = 1
  //
  // A = laplacian (eps) + mass (1/dt + exp (u))
  //     |______________| |_____________________|
  //           A0                     A1
  //
  // f = - mass ((1/dt) * (u - uold) + exp (u) - 1)  -  laplacian (eps) * u
  //     |_________________________________________|   |___________________|
  //                 f0                                         f1
  //
  // A * du = f
  // u = u + du

  std::vector<double> uold, uvold, uguess;
  uold.assign (n_nodes, 1.0);
  uvold.assign (n_nodes, 1.0);
  uguess.assign (n_nodes, 1.0);
  
  std::vector<double> u (uold);
  std::vector<double> du (uold);
  du.assign (n_nodes, 0.0);
  
  auto iu  = u.begin ();
  auto idu = u.begin ();
  auto iuo = uold.begin ();

  auto rhsfun = [&iu, &iuo, &dt] ()
    {return -(*iu - *(iuo++)) / dt - std::exp (*(iu++)) + 1.0; };

  double dtinv = 1 / dt;
  auto expudtinv = [&iu, &dtinv] ()
    { return std::exp (*(iu++)) + dtinv; };

  /// Allocate system matrix
  mumps *lin_solver = new mumps ();

  std::vector<double> xa;
  std::vector<int> ir, jc;
  
  sparse_matrix A;
  A.resize (n_nodes);
  dirichlet_bcs bcs;
  for (int i = 0; i < 4; ++i)
    bcs.push_back (std::make_tuple (0, i, [](double x, double y){return .0;}));
      
  bim2a_structure (tmsh, A);
  A.aij (xa, ir, jc, lin_solver->get_index_base ());

  std::vector<double> ecoeff (n_nodes);
  std::vector<double> ncoeff (n_elements);
  std::vector<double> f (n_nodes);

  ecoeff.assign (n_elements, 1.0);
  ncoeff.assign (n_nodes, 1.0);
  f.assign (n_nodes, 0.0);
 
  lin_solver->set_lhs_distributed ();
  lin_solver->set_distributed_lhs_structure (A.rows (), ir, jc);
  lin_solver->analyze ();
  
  int isave = 0;
  t_vect.push_back (t);
  t += dt;
  int sentinella = 0;
  for (auto t_save_p = t_save.begin (); t_save_p != t_save.end (); ++t_save_p)
    {

      while (told < (*t_save_p))
        {
          if (rank == 0)
            {
              if (t > (*t_save_p)) t = (*t_save_p);
	      dt = t - told;
	      
            }
	  if (rank == 0)
	    std::cout << "TIME : "<< t << std::endl;
	    
          MPI_Bcast (&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Bcast (&dtinv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Barrier (MPI_COMM_WORLD);

	  if (told != tvold)
	    for (int i = 0; i < n_nodes ; ++i)
	      uguess[i] = (t - tvold ) /
		(told - tvold) * uold[i] +
		dt / (tvold - told) * uvold[i];                              
	    
	  if (!all_non_negative (uguess))
	    {
	      if (rank == 0)
		std::cout << "Negative guess" <<std::endl;
	      for (int i = 0; i < n_nodes; ++i)
		uguess[i] = std::max (0.0, uguess[i]);
	    }

	  while (true) //////////////////////////////////////////////////////////////////
	    {
	      for (it_nonlin = 0; it_nonlin < MAX_IT; ++it_nonlin)
		{
		  
		  if (rank == 0) tic ();
		  f.assign (n_nodes, 0.0);
		  A.reset ();

		  // A0
		  ecoeff.assign (n_elements, eps_u);
		  ncoeff.assign (n_nodes, 1.0);
		  bim2a_advection_diffusion (tmsh, ecoeff, ncoeff, A);
              
		  // f1
		  sparse_matrix::col_iterator ja;
		  for (unsigned int ia = 0; ia < A.size (); ++ia)
		    if (A[ia].size ())
		      for (ja = A[ia].begin (); ja != A[ia].end (); ++ja)
			f[ia] -= A.col_val (ja) * u[A.col_idx (ja)];

		  // f0
		  ecoeff.assign (n_elements, 1.0);
		  iu = u.begin ();
		  iuo = uguess.begin ();
		  
		  std::generate (ncoeff.begin (), ncoeff.end (), rhsfun);
		  bim2a_rhs (tmsh, ecoeff, ncoeff, f);

		  // A1
		  ecoeff.assign (n_elements, 1.0);
		  iu = u.begin ();
		  iuo = uguess.begin ();
		  dtinv = 1 / dt;
	          std::generate (ncoeff.begin (), ncoeff.end (), expudtinv);
		  bim2a_reaction (tmsh, ecoeff, ncoeff, A);             
		  MPI_Barrier (MPI_COMM_WORLD);
		  if (rank == 0) toc ("assembly");
              
		  // TODO: MPI_Reduce f on rank 0

		  if (rank == 0) tic ();
		  std::copy (f.begin (), f.end (), du.begin ());
		  lin_solver->set_rhs (du);

		  // 	  bim2a_dirichlet_bc (tmsh, bcs, A, du);
        

		  A.aij_update (xa, ir, jc, lin_solver->get_index_base ());

		  lin_solver->set_distributed_lhs_data (xa);

		  lin_solver->factorize ();
		  lin_solver->solve ();
		  MPI_Barrier (MPI_COMM_WORLD);
		  if (rank == 0) toc ("solve");

		  if (rank == 0)
		    {
		      tic ();
		      residual_norm = 0.0;
		      std::for_each (du.begin (), du.end (), compute_norm);
		      residual_norm = std::sqrt (residual_norm);
		    }
		  MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              
		  for (int i = 0; i < n_nodes; ++i)
		    u[i] += du[i];
		  
		  MPI_Barrier (MPI_COMM_WORLD);
		  if (rank == 0) toc ("increment");
              
		  if (rank == 0)
		    std::cout << "iteration = "
			      << it_nonlin
			      << " incr norm = "
			      << residual_norm
			      << std::endl;
		  if (residual_norm <= MIN_RESIDUAL) break;
                        
		}
	      if (residual_norm <= MIN_RESIDUAL && all_non_negative (u))
	        break;
	      else
		{
		  if (rank == 0)
		    {
		      dt *= 0.5;
		      dtinv = 1 / dt;
		      t = told + dt;
		      std::cout << "Reducing time step : t = " << t << ", dt = "<< dt << std::endl;
		    }
		  MPI_Bcast (&dtinv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  MPI_Barrier (MPI_COMM_WORLD);

		  if (told != tvold)
		    {
		      for (int i = 0; i < n_nodes ; ++i)
			uguess[i] = (t - tvold ) /
			  (told - tvold) * uold[i] +
			  dt / (tvold - told) * uvold[i];                              
	    
		      if (!all_non_negative (uguess))
			{
			  if (rank == 0)
			    std::cout << "Negative guess" <<std::endl;
			  for (int i = 0; i < n_nodes; ++i)
			    uguess[i] = std::max (0.0, uguess[i]);
			    
			}
		    }
		  std::copy (uold.begin (), uold.end (), u.begin ());
         
		}
	    }
	  
          std::copy (uold.begin (), uold.end (), uvold.begin ());
          std::copy (u.begin (), u.end (), uold.begin ());
     	  tvold = told;
          told = t;
	  t_vect.push_back (t);
	  if (rank == 0)
	    {
	      if (residual_norm < 1e-14)
		dt *= 2;
	    
	      else
		dt = dt * std::min (std::sqrt(.38) * std::sqrt (MIN_RESIDUAL / residual_norm), 2.0);
	      
	      dt = std::min (dt ,  ((T - T0) / NT));
	      t += dt;
	      std::cout <<"----dt = "<< dt << std::endl;    
	      std::cout <<"----final step error = "<< residual_norm << std::endl;    
	    }
	  MPI_Bcast (&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Barrier (MPI_COMM_WORLD);

	}

      
      // std::cout << "SALVO AL TEMPO "<< t <<std::endl;

      //      if (rank == 0) tic ();

        tmsh.octbin_export ((std::string ("tumor_growth_u_")
                             + std::to_string (isave)).c_str (), u);

        tmsh.octbin_export ((std::string ("tumor_growth_f_")
                             + std::to_string (isave)).c_str (), f);
        
        tmsh.octbin_export ((std::string ("tumor_growth_du_")
                             + std::to_string (isave++)).c_str (), du);

      MPI_Barrier (MPI_COMM_WORLD);  
      //      if (rank == 0) toc ("export");

    }
  std::cout <<" t_vect : " << std::endl; 
  for (int i = 0 ; i < t_vect.size () ; ++i)
    std::cout<<t_vect[i]<<std::endl;

  std::cout <<" t_save : " << std::endl;
  for (int i = 0 ; i < t_save.size () ; ++i)
    std::cout<<t_save[i]<<std::endl;
 
  print_timing_report ();
  lin_solver->cleanup ();
  MPI_Finalize ();

  return 0;
}



