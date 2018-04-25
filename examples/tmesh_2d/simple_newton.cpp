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


constexpr int NUM_REFINEMENTS           = 5; //5;
constexpr double MIN_RESIDUAL           = 1.e-6;
constexpr double MIN_RESIDUAL_TIME_STEP = 1;
constexpr int NT                        = 10;
constexpr int MAX_IT                    = 50;
constexpr double DT                     = 1.0e-2;

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
all_non_negative (std::vector<double> vect, tmesh::idx_t &first_node,
		  int &last_node )
{
  for (int i = first_node; i < last_node; ++i)
    if (vect[i] < 0)
      return false;
  return true; 
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
  tmesh::idx_t num_local_nodes  = tmsh.num_owned_nodes ();
  tmesh::idx_t first_node  = tmsh.lnodes->global_offset;
  int last_node = first_node + num_local_nodes;  
 
  std::cout <<"RANK "<<rank<< " has : "<< num_local_nodes <<" number of nodes" <<std::endl;
  MPI_Barrier (MPI_COMM_WORLD);
  std::cout <<"RANK "<<rank<< " has as first node : "<< first_node <<std::endl;
  MPI_Barrier (MPI_COMM_WORLD);

  /// Time stepping variables
  double t           = T0;
  double told        = T0;
  double tvold       = T0;
  double dt          = DT;
  double dtold       = DT;
  int nt             = NT;

  std::array<double, NT+1> t_save {t};
  std::vector<double> t_vect;
  double dt_tsave = ((T - T0) / NT);
  
  auto its = [&t, &dt_tsave] ()
    {double out = t; t += dt_tsave ; return out; };
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

  std::vector<double> uold, uvold;
  uold.assign (n_nodes, 1.0);
  uvold.assign (n_nodes, 1.0);
  
  std::vector<double> u (uold);
  std::vector<double> u_local (n_nodes), f_global (n_nodes), du_global (n_nodes);
  std::vector<double> du;
  
  du.assign (n_nodes, 0.0);
  du_global.assign (n_nodes, 0.0);
  
  auto iu  = u.begin ();
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
  if (rank==0)
    t_vect.push_back (t);
  
  u_local.assign (n_nodes, 0.0);
  for (int i = first_node; i < last_node; ++i)
    u_local[i] = uold[i];
  
  int flag_while_tsave, flag_tvold, flag_neg, flag_neg_global;
  
  for (auto t_save_p = t_save.begin (); t_save_p != t_save.end (); ++t_save_p)
    {
      if (rank == 0)
	flag_while_tsave = (t < (*t_save_p));
      MPI_Bcast (&flag_while_tsave, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Barrier (MPI_COMM_WORLD);
      
      while (flag_while_tsave)
        {
	  if (rank == 0)
            {
	      t += dt;
              if (t > (*t_save_p)) t = (*t_save_p);
	      dt = t - told;
	      dtinv = 1 / dt;
	      std::cout << "TIME : "<< t << std::endl;
	    }
	  
          MPI_Bcast (&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Bcast (&dtinv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  //MPI_Barrier (MPI_COMM_WORLD);

	  if (rank == 0) flag_tvold = (told != tvold);
	  MPI_Bcast (&flag_tvold, 1, MPI_INT, 0, MPI_COMM_WORLD);
	  MPI_Barrier (MPI_COMM_WORLD);
    	  
	  if (flag_tvold)
	    for (int i = first_node; i < last_node; ++i)
	      u_local[i] = (t - tvold) /
		(dtold) * uold[i] + 
		dt / (-dtold) * uvold[i];                              
	    

	  flag_neg = ! all_non_negative (u_local, first_node, last_node);
	  flag_neg_global = 0;
	  MPI_Allreduce (&flag_neg, &flag_neg_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  MPI_Barrier (MPI_COMM_WORLD);
          if (flag_neg_global)
	    {
	      if (rank == 0)
		std::cout << "Negative guess" <<std::endl;
	      for (int i = first_node; i < last_node; ++i)
		u_local[i] = std::max (0.0, u_local[i]);
	    }
	  u.assign (n_nodes, 0.0);
	  MPI_Allreduce (&u_local[0], &u[0], n_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        	 
	  while (true) //////////////////////////////////////////////////////////////////
	    {
	      for (it_nonlin = 0; it_nonlin < MAX_IT; ++it_nonlin)
		{
		  // if (rank == 0) tic ();
		  f.assign (n_nodes, 0.0);
		  A.reset ();

		  // A0
		  ecoeff.assign (n_elements, eps_u);
		  ncoeff.assign (n_nodes, 0.1); 
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
		  iuo = uold.begin ();
		  
		  std::generate (ncoeff.begin (), ncoeff.end (), rhsfun);

		  bim2a_rhs (tmsh, ecoeff, ncoeff, f);
		  
		  // A1
		  ecoeff.assign (n_elements, 1.0);
		  iu = u.begin ();
		  iuo = uold.begin ();
		  // dtinv = 1 / dt;
	          std::generate (ncoeff.begin (), ncoeff.end (), expudtinv);
		  bim2a_reaction (tmsh, ecoeff, ncoeff, A);             
		  MPI_Barrier (MPI_COMM_WORLD);
		  //if (rank == 0) toc ("assembly");

		  bim2a_dirichlet_bc (tmsh, bcs, A, f);
		
		  if (rank == 0)
		    du_global.assign (n_nodes, 0);
		  MPI_Reduce (&f[0], &du_global[0], n_nodes, MPI_DOUBLE, MPI_SUM, 0,
			      MPI_COMM_WORLD);

		  if (rank == 0)
		    lin_solver->set_rhs (du_global);

		  /*	// QUI STAMPO SOLO PER VERIFICA	  
		  MPI_Barrier (MPI_COMM_WORLD);/////////////////////
		  if (rank == 0)
		    {
		      std::cout << "rank 0 :" << std::endl;
		      for (int i = 0; i < n_nodes; ++i)
			std::cout<< du[i]<<std::endl;
		    }
		  MPI_Barrier (MPI_COMM_WORLD);
		  if (rank == 1)
		    {
		      std::cout << "rank 1 :" << std::endl;
		      for (int i = 0; i < n_nodes; ++i)
			std::cout<< du[i]<<std::endl;
		    }
		  MPI_Barrier (MPI_COMM_WORLD);
		  if (rank == 2)
		    {
		      std::cout << "rank 2 :" << std::endl;
		      for (int i = 0; i < n_nodes; ++i)
			std::cout<< du[i]<<std::endl;
		    }		
		  MPI_Barrier (MPI_COMM_WORLD);
		  */

		  A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
		  lin_solver->set_distributed_lhs_data (xa);

		  lin_solver->factorize ();
		  lin_solver->solve ();
		  MPI_Barrier (MPI_COMM_WORLD);
        
		  // if (rank == 0) toc ("solve");

		  if (rank == 0)  // da parallelizzare
		    {
		      //   tic ();
		      residual_norm = 0.0;
		      std::for_each (du_global.begin (), du_global.end (), compute_norm);
		      residual_norm = std::sqrt (residual_norm);
		    }
		  MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  MPI_Bcast (&du_global[0], n_nodes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  for (int i = first_node; i < last_node; ++i)
		    u_local[i] += du_global[i];
		  
		  MPI_Barrier (MPI_COMM_WORLD);
		  u.assign (n_nodes, 0.0);
		  MPI_Allreduce (&u_local[0], &u[0], n_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 
		  // if (rank == 0) toc ("increment");
                 
		  if (rank == 0)
		    std::cout << "iteration = "
			      << it_nonlin
			      << " incr norm = "
			      << residual_norm
			      << std::endl;
		  
		  flag_neg = ! all_non_negative (u_local, first_node, last_node);
		  flag_neg_global = 0 ;
		  MPI_Allreduce (&flag_neg, &flag_neg_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	          MPI_Barrier (MPI_COMM_WORLD);
	          if (flag_neg_global)   it_nonlin = MAX_IT;
		  
		  if (residual_norm <= MIN_RESIDUAL)  break;
		  
        	}

	      if (residual_norm <= MIN_RESIDUAL && it_nonlin < MAX_IT)
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
		  MPI_Bcast (&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  // MPI_Barrier (MPI_COMM_WORLD);

		  if (rank == 0) flag_tvold = (told != tvold);
		  MPI_Bcast (&flag_tvold, 1, MPI_INT, 0, MPI_COMM_WORLD);
		  MPI_Barrier (MPI_COMM_WORLD);

		   if (flag_tvold)
		    {
		      for (int i = first_node; i < last_node; ++i)
			u_local[i] = (t - tvold) /
			  (dtold) * uold[i] + 
			  dt / (-dtold) * uvold[i];                              
	  	    
		      flag_neg = ! all_non_negative (u_local, first_node, last_node);
		      flag_neg_global = 0;
		      MPI_Allreduce (&flag_neg, &flag_neg_global, 1, MPI_INT, MPI_SUM,
				     MPI_COMM_WORLD);
		      MPI_Barrier (MPI_COMM_WORLD);
		      if (flag_neg_global)
			{
			  if (rank == 0)
			    std::cout << "Negative guess" <<std::endl;
			  for (int i = first_node; i < last_node; ++i)
			    u_local[i] = std::max (0.0, u_local[i]);
	        	}
		    }
		  else 
		    {
		      u_local.assign (n_nodes, 0.0); // per sicurezza, ma sarebbe da togliere
		      for (int i = first_node; i < last_node; ++i)
			u_local[i] = uold[i];
		    }
		   MPI_Allreduce (&u_local[0], &u[0], n_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 
		}
	    }
	  
	  for (int i = first_node; i < last_node; ++i) 
	    du[i] = u_local[i] - uold[i];

	  if (rank == 0)
	    du_global.assign (n_nodes, 0.0);
	  MPI_Reduce (&du[0], &du_global[0], n_nodes, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	 
	  std::copy (uold.begin (), uold.end (), uvold.begin ());
          std::copy (u.begin (), u.end (), uold.begin ());
          // if (rank == 0)
	  //  {
	  tvold = told;
	  told = t;
	  // }
	  // MPI_Bcast (&tvold, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  // MPI_Bcast (&told, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	  if (rank == 0)
	    t_vect.push_back (t);
	  dtold = dt;

	  if (rank == 0)
	    {
	      residual_norm = 0.0;
	      std::for_each (du_global.begin (), du_global.end (), compute_norm); // da rendere parallelo
	      residual_norm = std::sqrt (residual_norm);

	      if (residual_norm == 0) 
		dt *= 2;
	    
	      else
		dt = dt * std::min (std::sqrt (.38)
				    * std::sqrt (MIN_RESIDUAL_TIME_STEP / residual_norm), 2.0);
	      	        
	      dt = std::min (dt , dt_tsave);
	      std::cout <<"----dt = "<< dt << std::endl;    
	      std::cout <<"----final step error = "<< residual_norm << std::endl;
	      
	    }
	  MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  if (rank == 0)
	    flag_while_tsave = (t < (*t_save_p));
	  MPI_Bcast (&flag_while_tsave, 1, MPI_INT, 0, MPI_COMM_WORLD);
	  MPI_Barrier (MPI_COMM_WORLD);

	}

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

  if (rank == 0)
    {
      std::cout <<"--- t_vect --- " << std::endl; 
      for (int i = 0 ; i < t_vect.size () ; ++i)
	std::cout<<t_vect[i]<<std::endl;
    }
  // print_timing_report ();
  lin_solver->cleanup ();
  MPI_Finalize ();

  return 0;
}



