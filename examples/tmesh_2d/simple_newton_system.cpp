/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <cmath>
#include <algorithm>

#include <bim_timing.h>
#include <mumps_class.h>
#include <tmesh.h>
#include <quad_operators.h>


constexpr int NUM_REFINEMENTS           = 3; //5;
constexpr double MIN_RESIDUAL           = 1.e-6;
constexpr double MIN_RESIDUAL_TIME_STEP = 1;
constexpr int NT                        = 10;
constexpr int MAX_IT                    = 50;
constexpr double DT                     = 1.0e-2;

constexpr double T0    = 0.;
constexpr double T     = 10.;
constexpr double eps_u = 1;
constexpr double eps_v = 1;

constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {0., 0., 1., 0.,  1., 1., 0., 1.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {1, 2, 3, 4, 1};

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

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

  tmesh::idx_t num_global_nodes = tmsh.num_global_nodes ();
  tmesh::idx_t num_local_quadrants = tmsh.num_local_quadrants ();
  tmesh::idx_t num_owned_nodes  = tmsh.num_owned_nodes ();
  tmesh::idx_t global_offset  = tmsh.lnodes->global_offset; 
  tmesh::idx_t first_node  = global_offset;
  int last_node = first_node + num_owned_nodes; 
  int num_global_nodes2 = 2 * num_global_nodes;

  for (int i = 0; i < size; ++i)
    {
      if (rank == i)
	{
	  std::cout <<"RANK "<<rank<< " has : "<< num_local_quadrants <<" number of local quadrants"
		    <<std::endl;
	  std::cout <<"RANK "<<rank<< " has as first node : "<< global_offset <<std::endl;
	}
      MPI_Barrier (MPI_COMM_WORLD);
    }

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
  double              residual_norm, residual_norm_loc;
  int                 it_nonlin;
  auto compute_norm = [&residual_norm_loc] (double x)
    { residual_norm_loc += std::pow (x, 2); };

 /// PDE arrays system
  // (1/dt) * mass * (u - uold) + laplacian (eps_u) * u + mass * (exp (u - v)) = 1
  // (1/dt) * mass * (v - vold) + laplacian (eps_v) * v - mass * (exp (u - v)) = 0
  //
  // Auu = laplacian (eps_u) + mass (1/dt + exp (u - v))
  //      |_________________| |________________________|
  //             A0                     A1
  //
  // Auv =  mass (- exp (u - v))
  //       |________________________|
  //                 A1
  //
  // Avu =  mass (- exp (u - v))
  //      |________________________|
  //                 A1
  //
  // Avv = laplacian (eps_v) + mass (1/dt - exp (u - v))
  //      |_________________| |________________________|
  //             A0                     A1
  //
  //
  // fu = - mass ((1/dt) * (u - uold) + exp (u - v) - 1)  -  laplacian (eps_u) * u
  //        |_________________________________________|      |___________________|
  //                 f0                                            f1
  //
  // fv = - mass ((1/dt) * (v - vold) - exp (u - v))  -  laplacian (eps_v) * v
  //        |______________________________________|    |_____________________|
  //                 f0                                            f1
  //
  // A * [du dv] = [fu fv] ,  x = [u v]
  // [u v] = [u v] + [du dv], dx = [du dv]

  std::vector<double> xold, xvold;
  xold.assign (num_global_nodes2, 1.0);
  xvold.assign (num_global_nodes2, 1.0);
 
  std::vector<double> x (xold);
  std::vector<double> x_local (num_global_nodes2);
  std::vector<double> dx (num_global_nodes2), dx_global (num_global_nodes2);
 
  dx.assign (num_global_nodes2, 0.0);
  dx_global.assign (num_global_nodes2, 0.0);

  auto iu  = x.begin ();
  auto iuo = xold.begin ();
  auto iv  = x.begin () + num_global_nodes;
  auto ivo = xold.begin () + num_global_nodes;

  auto rhsfun_u = [&iu, &iuo, &iv, &dt] ()
    {return -(*iu - *(iuo++)) / dt - std::exp (*(iu++) - *(iv++)) + 1.0; };

  auto rhsfun_v = [&iu, &iv, &ivo ,&dt] ()
    {return -(*iv - *(ivo++)) / dt + std::exp (*(iu++) - *(iv++)); };

  double dtinv = 1 / dt;
  auto expudtinv_u = [&iu, &iv, &dtinv] ()
    { return  std::exp (*(iu++) - *(iv++)) + dtinv; };

  auto expudtinv_v = [&iu, &iv, &dtinv] ()
    { return - std::exp (*(iu++) - *(iv++)) + dtinv; };

  auto expuv = [&iu, &iv, &dtinv] ()
    { return - std::exp (*(iu++) - *(iv++)); };

  /// Allocate system matrix
  mumps *lin_solver = new mumps ();

  std::vector<double> xa, xa_tot;
  std::vector<int> ir, jc, ir_tot, jc_tot;
  
  sparse_matrix A;
  A.resize (num_global_nodes);
  
  dirichlet_bcs bcs;
  for (int i = 0; i < 4; ++i)
    bcs.push_back (std::make_tuple (0, i, [](double x, double y){return .0;}));
      
  bim2a_structure (tmsh, A);
  A.aij (xa, ir, jc, lin_solver->get_index_base ());
  {
    std::vector<int> ir_temp (ir.size ()), jc_temp (ir.size ());
    for (int i = 0 ; i < ir.size (); ++i)
      {
	ir_temp[i] = ir[i] + num_global_nodes;
	jc_temp[i] = jc[i] + num_global_nodes;
      }
    ir_tot.insert (ir_tot.end (), ir.begin (), ir.end ());
    ir_tot.insert (ir_tot.end (), ir.begin (), ir.end ());
    ir_tot.insert (ir_tot.end (), ir_temp.begin (), ir_temp.end ());
    ir_tot.insert (ir_tot.end (), ir_temp.begin (), ir_temp.end ());
 
    jc_tot.insert (jc_tot.end (), jc.begin (), jc.end ());
    jc_tot.insert (jc_tot.end (), jc_temp.begin (), jc_temp.end ());
    jc_tot.insert (jc_tot.end (), jc.begin (), jc.end ());
    jc_tot.insert (jc_tot.end (), jc_temp.begin (), jc_temp.end ());
  }

  xa_tot.resize (ir_tot.size ());

  lin_solver->set_lhs_distributed ();
  lin_solver->set_distributed_lhs_structure (2 * A.rows (), ir_tot, jc_tot);
  lin_solver->analyze ();

  std::vector<double> ecoeff (num_global_nodes);
  std::vector<double> ncoeff (num_local_quadrants);
  std::vector<double> f (num_global_nodes);

  ecoeff.assign (num_local_quadrants, 1.0);
  ncoeff.assign (num_global_nodes, 1.0);
  f.assign (num_global_nodes, 0.0);

  int isave = 0;
  if (rank==0)
    t_vect.push_back (t);


  auto iu_local_first = x_local.begin () + global_offset;
  auto iu_local_last =  x_local.begin () + global_offset + num_owned_nodes;
  std::vector<double>::iterator iu_local;

  auto iv_local_first = x_local.begin () + global_offset + num_global_nodes;
  auto iv_local_last =  x_local.begin () + global_offset + num_owned_nodes + num_global_nodes;
  std::vector<double>::iterator iv_local;

  auto iuo_first = xold.begin () + global_offset;
  auto iuvo_first = xvold.begin () + global_offset;
  std::vector<double>::iterator iuvo;

  auto ivo_first = xold.begin () + global_offset + num_global_nodes;
  auto ivvo_first = xvold.begin () + global_offset + num_global_nodes;
  std::vector<double>::iterator ivvo;

  auto idu_first = dx.begin () + global_offset;
  std::vector<double>::iterator idu;

  auto idv_first = dx.begin () + global_offset + num_global_nodes;
  std::vector<double>::iterator idv;
  
  auto idu_global_first = dx_global.begin () + global_offset;
  auto idu_global_last =  dx_global.begin () + global_offset + num_owned_nodes;
  std::vector<double>::iterator idu_global;

  auto idv_global_first = dx_global.begin () + global_offset + num_global_nodes;
  auto idv_global_last =  dx_global.begin () + global_offset + num_owned_nodes + num_global_nodes;
  std::vector<double>::iterator idv_global;

  x_local.assign (num_global_nodes2, 0.0);

  
  for (int i = 0; i < size; ++i)
    {
      if (rank == i)
	{
	  std:: cout << "RANK = "<< rank << std::endl;
	  std::cout << "ir size : "<<ir.size ()<< std::endl;
	  std::cout << "ir_tot size : "<<ir_tot.size ()<< std::endl;
    	}
      MPI_Barrier (MPI_COMM_WORLD);
    }


  int flag_while_tsave, flag_tvold, flag_neg_u, flag_neg_v, flag_neg, flag_neg_global;
  
  for (auto t_save_p = t_save.begin (); t_save_p != t_save.end (); ++t_save_p)
    {
      while (true)
        {
	  if (rank == 0)
	    flag_while_tsave = (t >= (*t_save_p));
	  MPI_Bcast (&flag_while_tsave, 1, MPI_INT, 0, MPI_COMM_WORLD);
	  if (flag_while_tsave) break; 	  

	  if (rank == 0)
	    {
	      if (t + dt > (*t_save_p))
		dt = (*t_save_p) - t;
	    }
	  while (true)
	    {
	      if (rank == 0)
		{
		  t += dt;
		  std::cout << "TIME : "<< t << std::endl;
		}
	      MPI_Bcast (&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	      MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	      dtinv = 1 / dt;

	      if (rank == 0) flag_tvold = (told != tvold);
	      MPI_Bcast (&flag_tvold, 1, MPI_INT, 0, MPI_COMM_WORLD);
	      iuo = iuo_first;
	      iuvo = iuvo_first;
	      ivo = ivo_first; 
	      ivvo = ivvo_first;
	      iv_local = iv_local_first;
	      if (flag_tvold)
		for ( iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
		  {
		  (*iu_local) = (t - tvold) /
		    (dtold) * (*(iuo++)) + 
		    dt / (-dtold) * (*(iuvo++));                              

		  (*(iv_local++)) = (t - tvold) /
		    (dtold) * (*(ivo++)) + 
		    dt / (-dtold) * (*(ivvo++));             
		  }
	      else
		for (iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
		  {
		    (*iu_local) = (*(iuo++));
		    (*(iv_local++)) = (*(ivo++));
		  }
	      flag_neg_u = any_of (iu_local_first, iu_local_last, [] (double ii) {return ii < 0;});
	      flag_neg_v = any_of (iv_local_first, iv_local_last, [] (double ii) {return ii < 0;});

	      if (flag_neg_u)
		for (iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
		  (*iu_local) = std::max (0.0, (*iu_local));
        	
	      if (flag_neg_v)
	        for (iv_local = iv_local_first; iv_local != iv_local_last; ++iv_local)
		  (*iv_local) = std::max (0.0, (*iv_local));
        	
	      MPI_Allreduce (&x_local[0], &x[0], num_global_nodes2, MPI_DOUBLE, MPI_SUM,
			     MPI_COMM_WORLD);
	      
	      for (it_nonlin = 0; it_nonlin < MAX_IT; ++it_nonlin)
		{

		  // Auu
		  
		  // if (rank == 0) tic ();
		  f.assign (num_global_nodes, 0.0);
		  A.reset ();

		  // A0
		  ecoeff.assign (num_local_quadrants, eps_u);
		  ncoeff.assign (num_global_nodes, 0.1); 
                  bim2a_advection_diffusion (tmsh, ecoeff, ncoeff, A);

		  // f1
		  sparse_matrix::col_iterator ja;
		  for (unsigned int ia = 0; ia < A.size (); ++ia)
		    if (A[ia].size ())
		      for (ja = A[ia].begin (); ja != A[ia].end (); ++ja)
			f[ia] -=  A.col_val (ja) * x[A.col_idx (ja)]; 

		  // f0
		  ecoeff.assign (num_local_quadrants, 1.0);
		  iu = x.begin ();
		  iuo = xold.begin ();
        	  iv = x.begin () + num_global_nodes ;
		  std::generate (ncoeff.begin (), ncoeff.end (), rhsfun_u);

		  bim2a_rhs (tmsh, ecoeff, ncoeff, f);
		  
		  // A1
		  ecoeff.assign (num_local_quadrants, 1.0);
		  iu = x.begin ();
		  iuo = xold.begin ();
		  iv = x.begin () + num_global_nodes ;
	          std::generate (ncoeff.begin (), ncoeff.end (), expudtinv_u);
		  bim2a_reaction (tmsh, ecoeff, ncoeff, A);             
		  // MPI_Barrier (MPI_COMM_WORLD);
		  //if (rank == 0) toc ("assembly");

		  // bim2a_dirichlet_bc (tmsh, bcs, A, f);

		  MPI_Reduce (&f[0], &dx_global[0], num_global_nodes, MPI_DOUBLE, MPI_SUM, 0,
			      MPI_COMM_WORLD);
		  
		  A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
		  for (int i = 0; i < xa.size (); ++i)
		    xa_tot[i] = xa[i];

		  // Auv and  Avu
		  
		  // if (rank == 0) tic ();
		  f.assign (num_global_nodes, 0.0);
		  A.reset ();
		  ncoeff.assign (num_global_nodes, 0.1); 

		  // A1
		  ecoeff.assign (num_local_quadrants, 1.0);
		  iu = x.begin ();
        	  iv = x.begin () + num_global_nodes ;
	         
		  std::generate (ncoeff.begin (), ncoeff.end (), expuv);
		  bim2a_reaction (tmsh, ecoeff, ncoeff, A);             
		  // MPI_Barrier (MPI_COMM_WORLD);
		  //if (rank == 0) toc ("assembly");

		  // bim2a_dirichlet_bc (tmsh, bcs, A, f);

		  
		  A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
		  for (int i = 0; i < xa.size (); ++i)
		    xa_tot[i + xa.size ()] = xa[i];
		  for (int i = 0; i < xa.size (); ++i)
		    xa_tot[i + 2 * xa.size ()] = xa[i];
		  // Avv
		  
		  // if (rank == 0) tic ();
		  f.assign (num_global_nodes, 0.0);
		  A.reset ();

		  // A0
		  ecoeff.assign (num_local_quadrants, eps_v);
		  ncoeff.assign (num_global_nodes, 0.1); 
                  bim2a_advection_diffusion (tmsh, ecoeff, ncoeff, A);

		  // f1
		  for (unsigned int ia = 0; ia < A.size (); ++ia)
		    if (A[ia].size ())
		      for (ja = A[ia].begin (); ja != A[ia].end (); ++ja)
			f[ia] -=  A.col_val (ja) * x[A.col_idx (ja) + num_global_nodes]; 
		  
		  // f0
		  ecoeff.assign (num_local_quadrants, 1.0);
		  iv = x.begin () + num_global_nodes;
		  ivo = xold.begin () + num_global_nodes;
        	  iu = x.begin ();
		  std::generate (ncoeff.begin (), ncoeff.end (), rhsfun_v);

		  bim2a_rhs (tmsh, ecoeff, ncoeff, f);
		  
		  // A1
		  ecoeff.assign (num_local_quadrants, 1.0);
		  iv = x.begin () + num_global_nodes;
		  ivo = xold.begin () + num_global_nodes;
        	  iu = x.begin ();
		  std::generate (ncoeff.begin (), ncoeff.end (), expudtinv_v);
		  bim2a_reaction (tmsh, ecoeff, ncoeff, A);             
		  // MPI_Barrier (MPI_COMM_WORLD);
		  //if (rank == 0) toc ("assembly");

		  //	  bim2a_dirichlet_bc (tmsh, bcs, A, f);

		  
		  MPI_Reduce (&f[0], &dx_global[num_global_nodes], num_global_nodes,
			      MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		  
		  A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
		  for (int i = 0; i < xa.size (); ++i)
		    xa_tot[i + 3 * xa.size ()] = xa[i];
		  
		  lin_solver->set_distributed_lhs_data (xa_tot);
		  if (rank == 0)
		    lin_solver->set_rhs (dx_global);

		  lin_solver->factorize ();
		  lin_solver->solve ();
		  // MPI_Barrier (MPI_COMM_WORLD);
       		  // if (rank == 0) toc ("solve");

		  residual_norm_loc = 0.0;
		  MPI_Bcast (&dx_global[0], num_global_nodes2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  /*
		  //*************
	      if (rank == 0)
		{
		  for (int i = 0; i<x.size (); ++i)
		    std::cout << dx_global[i]<<std::endl;
		  return 0;
		}
	      //************	      
	      */
		  std::for_each (idu_global_first, idu_global_last, compute_norm);
		  std::for_each (idv_global_first, idv_global_last, compute_norm);

		  MPI_Reduce (&residual_norm_loc, &residual_norm, 1, MPI_DOUBLE, MPI_SUM, 0,
				MPI_COMM_WORLD);
		  if (rank == 0)
		    residual_norm = std::sqrt (residual_norm);
		  MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  idu_global = idu_global_first;
		  for (iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
		    (*iu_local) += (*(idu_global++));

		  idv_global = idv_global_first;
		  for (iu_local = iv_local_first; iu_local != iv_local_last; ++iu_local)
		    (*iu_local) += (*(idv_global++));
		  
		  MPI_Allreduce (&x_local[0], &x[0], num_global_nodes2, MPI_DOUBLE, MPI_SUM,
				 MPI_COMM_WORLD);
 
		  // if (rank == 0) toc ("increment");
                 
		  if (rank == 0)
		    std::cout << "iteration = "
			      << it_nonlin
			      << " incr norm = "
			      << residual_norm
			      << std::endl;
		  
		  flag_neg = any_of (iu_local_first, iu_local_last, [] (double ii){return ii < 0;});
		  flag_neg += any_of (iv_local_first, iv_local_last, [] (double ii)
				      {return ii < 0;});

		  MPI_Allreduce (&flag_neg, &flag_neg_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
		      std::cout << "Reducing time step : t = " << t << ", dt = "<< dt << std::endl;
		    }
		}
	    }

	  idu = idu_first;
	  iuo = iuo_first;
	  for (iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
	    (*(idu++)) = (*(iu_local)) - (*(iuo++));
	  
	  idv = idv_first;
	  ivo = ivo_first;
	  for (iv_local = iv_local_first; iv_local != iv_local_last; ++iv_local)
	    (*(idv++)) = (*(iv_local)) - (*(ivo++));

	  MPI_Allreduce (&dx[0], &dx_global[0], num_global_nodes2, MPI_DOUBLE, MPI_SUM,
			 MPI_COMM_WORLD);
	 
	  std::copy (xold.begin (), xold.end (), xvold.begin ());
          std::copy (x.begin (), x.end (), xold.begin ());
          tvold = told;
	  told = t;
	 
	  if (rank == 0)
	    t_vect.push_back (t);
	  dtold = dt;

	  residual_norm_loc = 0.0;
	  std::for_each (idu_global_first, idu_global_last, compute_norm);
	  std::for_each (idv_global_first, idv_global_last, compute_norm);
	 
	  MPI_Reduce (&residual_norm_loc, &residual_norm, 1, MPI_DOUBLE, MPI_SUM, 0,
		      MPI_COMM_WORLD);
	  residual_norm = std::sqrt (residual_norm);
	     
	  if (rank == 0)
	    {
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

	}

      //      if (rank == 0) tic ();

      /*
        tmsh.octbin_export ((std::string ("tumor_growth_u_")
                             + std::to_string (isave)).c_str (), u);

        tmsh.octbin_export ((std::string ("tumor_growth_f_")
                             + std::to_string (isave)).c_str (), f);
        
        tmsh.octbin_export ((std::string ("tumor_growth_du_")
                             + std::to_string (isave++)).c_str (), du);
      */
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



