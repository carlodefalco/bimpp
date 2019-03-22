#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators_3d.h>

#include <simple_connectivity_3d.h>

#include <cassert>
#include <limits>

// uniform_refinement:
//    returns 1 ----> all quadrants are refined
static int
uniform_refinement (tmesh_3d::quadrant_iterator q)
{ return 1; }

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);

  using idx_t = tmesh_3d::idx_t;
  using q1_vec = q1_vec<std::vector<double>>;
  using gradient3 = gradient3<std::vector<double>>;
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  // number of initial uniform refinement steps
  constexpr unsigned unif_refine_steps = 3;
  // number of adaptive refinement steps
  constexpr unsigned adapt_refine_steps = 10;
  
  std::vector<idx_t>    nnodes;         // number of nodes
  std::vector<double>   h_step;         // mesh size

  // Problem parameters
  double epsilon = 1e-6;                // diffusion coefficient
  double theta = M_PI / 4;              // angle
  double n_coeff = 1/std::sqrt(3.0);    // normalization coefficient

  // Mesh generation
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  // Initial level of uniform refinement
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < unif_refine_steps; ++cycle)
    {
      tmsh.refine (recursive, partforcoarsen);
      tmsh.set_refine_marker (uniform_refinement);  
    }
  
  // Export mesh
  tmsh.vtk_export ("p4est_adr_test_1_metrics_initial_mesh");
  
  // Adaptive refinement loop
  for (int adapt = 0; adapt < adapt_refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " (rank " 
                << rank << ") ***" << std::endl;
      
      // Assemble matrix.
      sparse_matrix A;
      A.resize(tmsh.num_global_nodes());

      std::vector<double> alpha(tmsh.num_local_quadrants (), epsilon);
  
      q1_vec psi(tmsh.num_global_nodes (),0);
      //q1_vec psi(tmsh.num_owned_nodes ());
      //bim3a_solution_with_ghosts(tmsh,psi);
      double x = 0, y = 0, z = 0;
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int ii = 0; ii < 8; ++ii)
            {
              if (! quadrant->is_hanging (ii))
                {
                  x = quadrant->p(0, ii);
                  y = quadrant->p(1, ii);
                  z = quadrant->p(2, ii);
                  
                  psi[quadrant->gt(ii)] = (- x - y - z) / n_coeff;
                  psi[quadrant->gt(ii)] /= epsilon;
                }
            }
        }
      //psi.assemble(replace_op);

      // Reduce coefficients.

      q1_vec global_psi(tmsh.num_global_nodes(),0);
      MPI_Allreduce(psi.data(), global_psi.data(), 
                    psi.size(), MPI_DOUBLE, MPI_MAX, mpicomm);
			/*
      q1_vec global_psi(tmsh.num_owned_nodes ());
      bim3a_solution_with_ghosts(tmsh,global_psi);   
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int ii = 0; ii < 8; ++ii)
            if (! quadrant->is_hanging (ii))
              global_psi[quadrant->gt(ii)] = psi[quadrant->gt(ii)];
        }
      global_psi.assemble(replace_op);  
      */ 

      bim3a_advection_diffusion (tmsh, alpha, global_psi, A);
      
      // Assemble right-hand side.
      q1_vec rhs(tmsh.num_global_nodes(),0);
      //q1_vec rhs(tmsh.num_owned_nodes ());
      //bim3a_solution_with_ghosts(tmsh,rhs);  
      
      std::vector<double> f(tmsh.num_local_quadrants (), 0);

      q1_vec g(tmsh.num_global_nodes(),0);
      //q1_vec g(tmsh.num_owned_nodes ());
      //bim3a_solution_with_ghosts(tmsh,g);  
      
      bim3a_rhs (tmsh, f, g, rhs);

      // Set boundary conditions.
      func3 u0  = [] (double x, double y, double z) { return 0; };
      func3 u10 = [] (double x, double y, double z) 
      	{
      		double d = 1.;
      		if (x >= 0.5 && y >= (- x + 1.5))
      			d = 0.;
      		return d; 
      	};
      
      dirichlet_bcs3 bcs;
      bcs.push_back (std::make_tuple(0, 0, u0));
      bcs.push_back (std::make_tuple(0, 1, u0));
      bcs.push_back (std::make_tuple(0, 2, u0));
      bcs.push_back (std::make_tuple(0, 3, u0));
      bcs.push_back (std::make_tuple(0, 4, u0));
      bcs.push_back (std::make_tuple(0, 5, u10));
      
      bim3a_dirichlet_bc (tmsh, bcs, A, rhs);
      
      // Solve problem.
      std::cout << "Solving linear system (rank "
                << rank << ")" << std::endl;

      mumps mumps_solver(true);
      
      std::vector<double> vals;
      std::vector<int> irow, jcol;
      
      A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
      mumps_solver.set_distributed_lhs_data (vals);
      
      // Reduce rhs (so that rank 0 has the actual rhs).
      
      q1_vec global_rhs(tmsh.num_global_nodes(),0);
      MPI_Reduce(rhs.data(), global_rhs.data(), 
                  rhs.size(),MPI_DOUBLE, MPI_SUM, 0, mpicomm);
			
      /*
      q1_vec global_rhs(tmsh.num_owned_nodes ());
      bim3a_solution_with_ghosts(tmsh,global_rhs);   
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int ii = 0; ii < 8; ++ii)
            if (! quadrant->is_hanging (ii))
              global_rhs[quadrant->gt(ii)] = rhs[quadrant->gt(ii)];
        }
      global_rhs.assemble(); 
      */

      MPI_Barrier(mpicomm);
      if (rank == 0)
        mumps_solver.set_rhs(global_rhs);

      // Solve.
      std::cout << "\tanalyze (rank " << rank << ")" << std::endl;
      int analyze_res =  mumps_solver.analyze ();
      std::cout << "\tanalyze_res = " << analyze_res 
      					<< " (rank " << rank << ")" << std::endl;
      std::cout << "\tfactorize (rank " << rank << ")" << std::endl;
			mumps_solver.factorize ();
      std::cout << "\tsolve (rank " << rank << ")" << std::endl;
			mumps_solver.solve ();
      std::cout << "\tcleanup (rank " << rank << ")" << std::endl;
      mumps_solver.cleanup ();

      // Export solution.

      MPI_Bcast(global_rhs.data(), global_rhs.size(), MPI_DOUBLE, 0, mpicomm);

      tmsh.octbin_export ((std::string("p4est_adr_test_1_metrics_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::cout << "Done (rank " << rank << ")" << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient and estimator (rank "
      					<< rank << ")" << std::endl;
      
      std::cout << "\tgradient (rank " << rank << ")" << std::endl;
      gradient3 du = bim3c_quadtree_pde_recovered_gradient(tmsh, global_rhs);
      std::cout << "\tsolution (rank " << rank << ")" << std::endl;
      q2_vec3 u_star = bim3c_quadtree_pde_recovered_solution(tmsh, global_rhs, 
                                                              du);
      
      // Export reconstructed gradient
      tmsh.octbin_export ((std::string("p4est_adr_test_1_metrics_du_x_")
                           + std::to_string(adapt)).c_str(), std::get<0>(du));
      tmsh.octbin_export ((std::string("p4est_adr_test_1_metrics_du_y_")
                           + std::to_string(adapt)).c_str(), std::get<1>(du));
      tmsh.octbin_export ((std::string("p4est_adr_test_1_metrics_du_z_")
                           + std::to_string(adapt)).c_str(), std::get<2>(du));
      
      // Solution estimator
      auto estimator = [& u_star, & global_rhs] (tmesh_3d::quadrant_iterator q)
        { return estimator_sol (q, u_star, global_rhs); };
      
      std::cout << "\tmetrics (rank " << rank << ")" << std::endl;
      
      // Set marker for refinement
      double tol = 1e-3;
      tmsh.set_metrics_marker (estimator, tol, 4);
      
      // Compute metrics and h.
      std::vector<double> metrics(tmsh.num_local_quadrants ());
      
      double  hx = 0, hy = 0, hz = 0,
              h = std::numeric_limits<double>::max (),
              global_h = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          metrics[quadrant->get_forest_quad_idx ()] =
            estimator (quadrant) * std::sqrt (tmsh.num_global_quadrants ())
            / tol;
          
          hx = quadrant->p(0, 7) - quadrant->p(0, 0);
          hy = quadrant->p(1, 7) - quadrant->p(1, 0);
          hz = quadrant->p(2, 7) - quadrant->p(2, 0);
          
          h = std::min(h, std::sqrt(hx*hx + hy*hy + hz*hz));
        }
      
      // Export metrics
      tmsh.octbin_export_quadrant ((std::string("p4est_adr_test_1_metrics_hx_")
                                    + std::to_string(adapt)).c_str(), metrics);
      
      // Compute global mesh size
      MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      h_step.push_back (global_h);
      
      std::cout << "Done (rank " << rank << ")\n" << std::endl;
      
      // Break if the number of global nodes is too large
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Refine.
      tmsh.metrics_refine (1e3);
      
      // Export new mesh
      tmsh.vtk_export ((std::string("p4est_adr_test_1_metrics_newmesh_")
                        + std::to_string(adapt)).c_str());
    }
  
  if (rank == 0)
    for (unsigned step = 0; step < nnodes.size(); ++step)
      std::cout << "Step " << step << ", #nodes: "
                << nnodes[step] << ", h: "
                << h_step[step] << std::endl;
  
  MPI_Finalize ();
  
  return 0;
}
