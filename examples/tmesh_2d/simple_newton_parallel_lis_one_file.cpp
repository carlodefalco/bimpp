/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <cmath>
#include <bim_timing.h>
#include <mumps_class.h>
#include <lis_class.h>
#include <lis_distributed_class.h>
#include <linear_solver.h>    
#include <tmesh.h>
#include <quad_operators.h>
#include <lis.h>
#include <string>
#include <bim_sparse.h>

#include <bim_sparse_distributed.h>
#include <algorithm>
#include <sstream>


constexpr int NUM_REFINEMENTS           = 2; //5;
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

struct
non_local_t
{
  std::vector<int> row_ptr, col_ind;
  std::vector<double> a;
  
  void
  csr (distributed_sparse_matrix &A, int is, int ie);
  // TA : ho aggiunto questo metodo che aggiorna i valori di a,
  // dato che row_ptr e col_ind rimangono costanti
  void
  csr_update (distributed_sparse_matrix &A, int is, int ie);
};

void
non_local_t::csr (distributed_sparse_matrix &A, int is, int ie)
{
  A.set_properties ();
  a.reserve (A.nnz); col_ind.reserve (A.nnz);
  row_ptr.resize (A.rows () + 1);

  int idx = 0, idr = 0;

  distributed_sparse_matrix::col_iterator jj;
  
  for (size_t ii = 0; ii < A.size (); ++ii)
    {
      row_ptr[idr++] = idx;
      if ((A[ii].size () > 0) && ((ii < is) || (ii >= ie)))
        {
          for (jj  = A[ii].begin (); jj != A[ii].end (); ++jj)
            {
              col_ind.push_back (A.col_idx (jj));
              a.push_back (A.col_val (jj));
              idx++;
            }
        }
    }
  std::fill (row_ptr.begin () + idr, row_ptr.end (), col_ind.size ());
};  
void
non_local_t::csr_update (distributed_sparse_matrix &A, int is, int ie)
{
  int idx = 0;
  for (int i = 1; i < row_ptr.size (); ++i)
    for (int j = 0; j < row_ptr[i] - row_ptr[i-1]; ++j)
      {
	a[idx] = A[i-1][col_ind[idx]];
	idx ++;
      }
};  

non_local_t non_local;

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
  uold.assign (num_global_nodes, 1.0);
  uvold.assign (num_owned_nodes, 1.0);

  // uold e u hanno dimensione piena, cioè il numero totale di nodi,
  // perchè servono intere per il calcolo di A e f con bim2a
  std::vector<double> u (num_global_nodes, 1.0);
  std::vector<double> u_local (num_owned_nodes);
  std::vector<double>  du_local (num_owned_nodes);
 
  du_local.assign (num_owned_nodes, 0.0);

  auto iu  = u.begin ();
  auto iuo = uold.begin ();
  auto rhsfun = [&iu, &iuo, &dt] ()
    {return -(*iu - *(iuo++)) / dt - std::exp (*(iu++)) + 1.0; };

  double dtinv = 1 / dt;
  auto expudtinv = [&iu, &dtinv] ()
    { return std::exp (*(iu++)) + dtinv; };

  std::vector<double> xa;
  std::vector<int> ir, jc;
  distributed_sparse_matrix A;
  A.resize (num_global_nodes);
  dirichlet_bcs bcs;
  for (int i = 0; i < 4; ++i)
    bcs.push_back (std::make_tuple (0, i, [](double x, double y){return .0;}));
      
  bim2a_structure (tmsh, A);
 
  std::vector<double> ecoeff (num_global_nodes);
  std::vector<double> ncoeff (num_local_quadrants);
  std::vector<double> f (num_global_nodes);

  ecoeff.assign (num_local_quadrants, 1.0);
  ncoeff.assign (num_global_nodes, 1.0);
  f.assign (num_global_nodes, 0.0);
  
  /// #################################################################################
  ///  TA: inizia la parte da mettere nella classe lis_distributed
  LIS_SOLVER solver;
  LIS_MATRIX A_lis;
  LIS_VECTOR b, x_lis;
  LIS_INT iter;
  LIS_INT *row,  *col;
  LIS_SCALAR *value;
  double time;
  bool have_initial_guess = false, initialized = false;
  int n, n_row;
  int max_iter = 1000;
  int restart_iterations = 40;
  double tolerance = 1.0e-12;
  std::string  iterative_method = "bicg";
  std::string  preconditioner = "none";
  std::string  convergence_condition = "nrm2_r";
  std::string  options_iterative_method = " ";
  std::string  option_string = "";
  bool  option_string_set = false ;
  bool  verbose = true;
  static const int index_base = 0;
  double *rhs, *data, *initial_guess;
  std::vector<int> map_n (size), map_row_s (size);
	
  n_row =  num_owned_nodes;
  // TA: creo solver, matrice e vettori all'inizio e poi li distruggo solo alla fine di tutto	  
  lis_matrix_create (MPI_COMM_WORLD, &A_lis);
  lis_vector_create (MPI_COMM_WORLD, &b);
  lis_vector_create (MPI_COMM_WORLD, &x_lis);
  lis_solver_create (&solver);
  lis_matrix_set_size (A_lis, num_owned_nodes, 0);
  LIS_INT is, ie;
  lis_matrix_get_range (A_lis, &is, & ie);
  int nloc = ie - is;
  MPI_Allgather(&is, 1, MPI_INT, &(map_row_s[0]), 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather (&n_row, 1, MPI_INT, &(map_n[0]), 1, MPI_INT, MPI_COMM_WORLD);
 
  A.set_ranges (is, ie);
  //A.remap ();
  // A.assemble ();
  lis_vector_set_size (b,  num_owned_nodes, 0);
  lis_vector_duplicate (b, &x_lis);
  // ##########################################################################################
  int isave = 0;
  if (rank==0)
    t_vect.push_back (t);

  auto iu_local_first = u_local.begin ();
  auto iu_local_last =  u_local.begin () + num_owned_nodes;
  std::vector<double>::iterator iu_local;

  auto iuo_first = uold.begin () + global_offset;
  auto iuo_last = uold.begin () + global_offset + num_owned_nodes;
  
  auto iuvo_first = uvold.begin ();
  std::vector<double>::iterator iuvo;
  
  auto idu_local_first = du_local.begin ();
  auto idu_local_last =  du_local.begin () + num_owned_nodes;
  std::vector<double>::iterator idu_local;

  u_local.assign (num_owned_nodes, 0.0);

  int flag_while_tsave, flag_tvold, flag_neg, flag_neg_global;
  int count = 0;
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
	      if (flag_tvold)
		for ( iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
		  (*iu_local) = (t - tvold) /
		    (dtold) * (*(iuo++)) + 
		    dt / (-dtold) * (*(iuvo++));                              
	      else
		for (iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
		  (*iu_local) = (*(iuo++));

	      flag_neg = any_of (iu_local_first, iu_local_last, [] (double ii) {return ii < 0;});
	      if (flag_neg)
		{
		  if (rank == 0)
		    std::cout << "Negative guess" <<std::endl;
		  for (iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
		    (*iu_local) = std::max (0.0, (*iu_local));
		}
	      MPI_Allgatherv (&u_local[0], num_owned_nodes, MPI_DOUBLE, &u[0],
				 &map_n[0], &map_row_s[0], MPI_DOUBLE, MPI_COMM_WORLD);
        	 
	      for (it_nonlin = 0; it_nonlin < MAX_IT; ++it_nonlin)
		{

		  count ++;
		  // if (rank == 0) tic ();
		  f.assign (num_global_nodes, 0.0);
		  A.reset ();

		  // A0
		  ecoeff.assign (num_local_quadrants, eps_u);
		  ncoeff.assign (num_global_nodes, 0.1); 
                  bim2a_advection_diffusion (tmsh, ecoeff, ncoeff, A);
	 
		  // f1
		  distributed_sparse_matrix::col_iterator ja;
		  for (unsigned int ia = 0; ia < A.size (); ++ia)
		    if (A[ia].size ())
		      for (ja = A[ia].begin (); ja != A[ia].end (); ++ja)
			f[ia] -= A.col_val (ja) * u[A.col_idx (ja)]; 

	        
		  // f0
		  ecoeff.assign (num_local_quadrants, 1.0);
		  iu = u.begin ();
		  iuo = uold.begin ();
		  
		  std::generate (ncoeff.begin (), ncoeff.end (), rhsfun);

		  bim2a_rhs (tmsh, ecoeff, ncoeff, f);
		  
		  // A1
		  ecoeff.assign (num_local_quadrants, 1.0);
		  iu = u.begin ();
		  iuo = uold.begin ();
	          std::generate (ncoeff.begin (), ncoeff.end (), expudtinv);
		  bim2a_reaction (tmsh, ecoeff, ncoeff, A);             
		  // MPI_Barrier (MPI_COMM_WORLD);
		  //if (rank == 0) toc ("assembly");
	      
		  bim2a_dirichlet_bc (tmsh, bcs, A, f);
		 
		  for (int i = 0; i < size; ++i)
		    MPI_Reduce (&f[map_row_s[i]], &du_local[0], map_n[i], MPI_DOUBLE, MPI_SUM, i,
				MPI_COMM_WORLD);


		  if (count < 4)
		    {
		      A.remap ();
		      A.assemble ();
		      A.csr (xa, jc, ir, 0);
		    }
		  else
		    {
                      A.assemble ();
		      A.csr (xa, jc, ir, 0);
		      
		    }
	        
		  lis_matrix_create (MPI_COMM_WORLD, &A_lis);
		  lis_matrix_set_size (A_lis, num_owned_nodes, 0);

	          lis_matrix_set_csr (jc.size () , &ir[0], &jc[0], &xa[0], A_lis);
		  lis_matrix_assemble (A_lis);
		  MPI_Bcast (&have_initial_guess, 1, MPI_INT, 0, MPI_COMM_WORLD);
		  		      
		  
	          for (int i =0; i < num_owned_nodes; ++i)
		    {
		      lis_vector_set_value (LIS_INS_VALUE, i + is, du_local[i], b);
		      if (have_initial_guess)
			lis_vector_set_value (LIS_INS_VALUE, i + is,
					      initial_guess[i], x_lis);
		    }
        
		  if (size == 1)
		    lis_output_vector (b, LIS_FMT_MM, "b_giusta.mm");

		  if (size == 2)
		    lis_output_vector (b, LIS_FMT_MM, "b.mm");

		  if (size == 1)
		    lis_output_matrix (A_lis, LIS_FMT_MM, "A_giusta.mm");

		  if (size == 2)
		    lis_output_matrix (A_lis, LIS_FMT_MM, "A.mm");
		  
		  // SOLVE !
		  char* options = 0;
		  if (! option_string_set)
		    {
		      std::stringstream opt;
		      opt << "-maxiter " << max_iter
			  << " -restrart " << restart_iterations
			  << " -tol " << tolerance
			  << " -i " << iterative_method
			  << options_iterative_method
			  << " -p " << preconditioner
			  << " -conv_cond " << convergence_condition;
          
		      if (have_initial_guess)
			opt << " -initx_zeros false ";
		      else
			opt << " -initx_zeros true ";
		      
		      option_string = opt.str ();
		      option_string_set = true;
		    }
		  options = new char[option_string.length () + 1];
		  std::copy (option_string.begin (),
			     option_string.end (), options);
		  // std::cout << options << std::endl;
		  lis_solver_set_option (options, solver);

		  lis_solve (A_lis, b, x_lis, solver);
		  
		  lis_solver_get_iter (solver, &iter);
		  lis_solver_get_time (solver, &time);
		  
		  delete [] options;
		  
		  //gather solution vector
		  double temp = 0.0;
		  for (int i = 0 ; i < num_owned_nodes; ++i)
		    {
		      lis_vector_get_value (x_lis, i + is, &temp);
		      du_local[i] = temp;
		    }
        
		  if (size == 1)
		    lis_output_vector (x_lis, LIS_FMT_MM, "x_giusta.mm");
		  if (size == 2)
		    lis_output_vector (x_lis, LIS_FMT_MM, "x.mm");
		  
		  		  
		  if (verbose && rank == 0)
		    std::cout << std::endl
			      << "Number of iterations = " << iter
			      << std::endl
			      << "Elapsed time = " << time << std::endl;

		  MPI_Barrier (MPI_COMM_WORLD);
       		  // if (rank == 0) toc ("solve");
		  
		  residual_norm_loc = 0.0;
		  std::for_each (idu_local_first, idu_local_last, compute_norm);
		  MPI_Reduce (&residual_norm_loc, &residual_norm, 1, MPI_DOUBLE, MPI_SUM, 0,
			      MPI_COMM_WORLD);
		  if (rank == 0)
		    residual_norm = std::sqrt (residual_norm);
		  MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		  idu_local = idu_local_first;
		  for (iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
		    (*iu_local) += (*(idu_local++));
		  

		  MPI_Allgatherv (&u_local[0], num_owned_nodes, MPI_DOUBLE, &u[0],
				  &map_n[0], &map_row_s[0], MPI_DOUBLE, MPI_COMM_WORLD);
        
		  // if (rank == 0) toc ("increment");
                 
		  if (rank == 0)
		    std::cout << "iteration = "
			      << it_nonlin
			      << " incr norm = "
			      << residual_norm
			      << std::endl;
		  
		
		  flag_neg = any_of (iu_local_first, iu_local_last, [] (double ii){return ii < 0;});
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

	  idu_local = idu_local_first;
	  iuo = iuo_first;
	  for (iu_local = iu_local_first; iu_local != iu_local_last; ++iu_local)
	    (*(idu_local++)) = (*(iu_local)) - (*(iuo++));
	 
	  std::copy (iuo_first , iuo_last, uvold.begin ());
          std::copy (u.begin (), u.end (), uold.begin ());
          tvold = told;
	  told = t;
	 
	  if (rank == 0)
	    t_vect.push_back (t);
	  dtold = dt;

	  residual_norm_loc = 0.0;
	  std::for_each (idu_local_first, idu_local_last, compute_norm);
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
                             + std::to_string (isave++)).c_str (), u);
      */
			     /*
        tmsh.octbin_export ((std::string ("tumor_growth_f_")
                             + std::to_string (isave)).c_str (), f);
        
        tmsh.octbin_export ((std::string ("tumor_growth_du_")
                             + std::to_string (isave++)).c_str (), du_local);
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
  //lin_solver->cleanup ();
  
  //lis_matrix_destroy (A_lis);
  //x_lis_solver_destroy (solver);
  //lis_vector_destroy (b);
  //  lis_vector_destroy (x_lis);
  MPI_Finalize ();

  return 0;
}
