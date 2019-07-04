#include <bim_sparse_distributed.h>
#include <mumps_class.h>
#include <quad_operators_3d.h>
#include <bim_timing.h>

#include <simple_connectivity_3d.h>

#include <cassert>
#include <limits>

// uniform_refinement:
//    returns 1 ----> all quadrants are refined
static int
uniform_refinement (tmesh_3d::quadrant_iterator q)
{ return 1; }

// Number of refinement steps
constexpr unsigned unif_refine_steps = 2;   // initial uniform refinement
constexpr unsigned adapt_refine_steps = 10;  // adaptive refinement

// Tolerance for refinement
constexpr double tol = 1e-3;

// main:
//
int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);

  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;

  using q1_vec = q1_vec<distributed_vector>;
  using gradient3 = gradient3<q1_vec>;
  using idx_t = tmesh_3d::idx_t;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  // Problem parameters
  constexpr double epsilon = 1e-10;     // diffusion coefficient
  double n_coeff = 1/std::sqrt(3.0);    // normalization coefficient
  
  // Mesh parameters
  std::vector<idx_t>    nnodes (adapt_refine_steps, 0);     // number of nodes
  std::vector<double>   h_step (adapt_refine_steps, 0.);    // mesh size

  // Estimators at every step
  std::vector<double> estSol (adapt_refine_steps,0.);  // ||u^* - u||_L^2(q)
  std::vector<double> estGrad (adapt_refine_steps,0.); // ||grad^* u - grad u||_L^2(q)

  // Mesh generation
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  // Initial uniform refinement
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < unif_refine_steps; ++cycle)
    {
      tmsh.refine (recursive, partforcoarsen);
      tmsh.set_refine_marker (uniform_refinement);
    }
  
  // Export initial mesh
  tmsh.vtk_export ("p4est_adr_test_1_metrics_initial_mesh");
  
  // Adaptive refinement loop
  for (int adapt = 0; adapt < adapt_refine_steps; ++adapt)
    {
      MPI_Barrier (mpicomm); 
      if (rank == 0)
        std::cout << "*** Step "  << adapt << " ***" << std::endl;
      
      // Compute coefficients

      // diffusion
      std::vector<double> alpha(tmsh.num_local_quadrants (), epsilon);
      q1_vec psi(tmsh.num_owned_nodes ());

      // rhs
      std::vector<double> f(tmsh.num_local_quadrants (), 0);
      q1_vec g(tmsh.num_owned_nodes ());

      double x = .0, y = .0, z = .0;
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

                  psi[quadrant->gt(ii)] = (x + y - z) * n_coeff / epsilon;
                  g[quadrant->gt(ii)] = 0.;
                }
              else
                for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                  {
                    psi[quadrant->gparent (jj, ii)] += 0.;
                    g[quadrant->gparent (jj, ii)] += 0.;
                  }
            }
        }
      psi.assemble(replace_op);
      g.assemble(replace_op);

      // Assemble matrix.
      distributed_sparse_matrix A;
      A.set_ranges(tmsh.num_owned_nodes());

      // advection_diffusion
      bim3a_advection_diffusion (tmsh, alpha, psi, A);

      // Assemble right-hand side.
      q1_vec rhs(tmsh.num_owned_nodes ());
      bim3a_rhs (tmsh, f, g, rhs);

      // Set boundary conditions.
      func3 u0  = [] (double x, double y, double z) { return 0; };
      func3 u10 = [] (double x, double y, double z) 
      {
        return ((x <= 0.5 && y <= (- x + 0.5)) ? 0. : 1.);
      };
      //
      dirichlet_bcs3 bcs;
      bcs.push_back (std::make_tuple(0, 0, u0));
      bcs.push_back (std::make_tuple(0, 1, u0));
      bcs.push_back (std::make_tuple(0, 2, u0));
      bcs.push_back (std::make_tuple(0, 3, u0));
      bcs.push_back (std::make_tuple(0, 4, u0));
      bcs.push_back (std::make_tuple(0, 5, u10));
      //
      bim3a_dirichlet_bc (tmsh, bcs, A, rhs);

      A.assemble();
      rhs.assemble();
      
      // Solve problem.
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }

      // Initialize MUMPS solver
      mumps mumps_solver;

      // Set distributed structure of rhs
      mumps_solver.set_rhs_distributed(rhs);
      
      // Set distributed structure of lhs
      std::vector<double> vals;
      std::vector<int> irow, jcol;
      
      A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);

      MPI_Barrier (mpicomm); if (rank == 0) { toc ("prepare solver "); }

      // Analyze
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      mumps_solver.analyze ();
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("analyze "); }

      // Set distributed data of lhs
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      mumps_solver.set_distributed_lhs_data (vals);
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("set lhs data "); }

      // Factorize
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
			mumps_solver.factorize ();
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("factorize "); }

      // Solve
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
			mumps_solver.solve ();
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("solve "); }

      // Get solution of linear system
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      q1_vec mumps_result = mumps_solver.get_distributed_solution();

      q1_vec result (tmsh.num_owned_nodes(), mpicomm);

      for (auto ii = result.get_range_start (); 
                ii != result.get_range_end (); 
                ++ii)
        result[ii] = mumps_result[ii];

      bim3a_solution_with_ghosts (tmsh, result, replace_op);
      result.assemble (replace_op);
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("scatter solution "); }

      // Export solution
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      tmsh.octbin_export ((std::string("p4est_adr_test_1_metrics_u_")
                           + std::to_string(adapt)).c_str(), result);
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("export solution "); }

      // Cleanup
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      mumps_solver.cleanup ();
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("cleanup "); }
      
      // Compute reconstructed gradient
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      gradient3 du = bim3c_quadtree_pde_recovered_gradient(tmsh, result);
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("gradient recovery "); }

      // Export reconstructed gradient
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      // 
      // derivative along x
      tmsh.octbin_export ((std::string ("p4est_adr_test_1_metrics_dudx_")
                          + std::to_string(adapt)).c_str(), std::get<0>(du));
      // 
      // derivative along y
      tmsh.octbin_export ((std::string ("p4est_adr_test_1_metrics_dudy_")
                          + std::to_string(adapt)).c_str(), std::get<1>(du));
      // 
      // derivative along z
      tmsh.octbin_export ((std::string ("p4est_adr_test_1_metrics_dudz_")
                          + std::to_string(adapt)).c_str(), std::get<2>(du));
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("export gradient "); }

      // Compute reconstructed solution
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      q2_vec3 u_star = bim3c_quadtree_pde_recovered_solution(tmsh, result, du);
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("solution recovery "); }
      
      // Define estimators
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }

      // solution estimator
      auto estimator = [& u_star, & result] (tmesh_3d::quadrant_iterator q)
        { return estimator_sol (q, u_star, result); };

      // gradient estimator
      auto grad_estimator = [& du, & result] (tmesh_3d::quadrant_iterator q)
        { return estimator_grad (q, du, result); };

      MPI_Barrier (mpicomm); if (rank == 0) { toc ("define estimators "); }
      
      // Compute metrics, mesh size and estimators
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }
      std::vector<double> metrics(tmsh.num_local_quadrants ());

      std::vector<double> est_sol (tmsh.num_local_quadrants ());
      std::vector<double> est_grad (tmsh.num_local_quadrants ());
      
      double  hx = 0, hy = 0, hz = 0,
              h = std::numeric_limits<double>::max ();

      double estsol = 0.0, estgrad = 0.0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          metrics[quadrant->get_forest_quad_idx ()] =
            (estimator (quadrant) * std::sqrt (tmsh.num_global_quadrants ()) 
            / tol);

          // ||u^* - u||_L^2(q)
          est_sol[quadrant->get_forest_quad_idx ()] = estimator(quadrant);
          estsol += std::pow(est_sol[quadrant->get_forest_quad_idx ()], 2);

          // ||grad^* u - grad u||_L^2(q)          
          est_grad[quadrant->get_forest_quad_idx ()] = grad_estimator(quadrant);
          estgrad += std::pow(est_grad[quadrant->get_forest_quad_idx ()], 2);
          
          hx = quadrant->p(0, 7) - quadrant->p(0, 0);
          hy = quadrant->p(1, 7) - quadrant->p(1, 0);
          hz = quadrant->p(2, 7) - quadrant->p(2, 0);
          
          h = std::min(h, std::sqrt(hx*hx + hy*hy + hz*hz));
        }
      
      // Export metrics
      tmsh.octbin_export_quadrant ((std::string("p4est_adr_test_1_metrics_hx_")
                                    + std::to_string(adapt)).c_str(), metrics);

      // Export estimators
      tmsh.octbin_export_quadrant ((std::string("p4est_adr_test_1_metrics_est_sol_")
                                    + std::to_string(adapt)).c_str(), est_sol);
      tmsh.octbin_export_quadrant ((std::string("p4est_adr_test_1_metrics_est_grad_")
                                    + std::to_string(adapt)).c_str(), est_grad);
      
      // Compute global mesh size
      MPI_Reduce(&h, &h_step[adapt], 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);

      // Compute global errors
      //
      // ||u^* - u||_L^2(q)
      MPI_Reduce (&estsol, &estSol[adapt], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      estSol[adapt] = std::sqrt(estSol[adapt]);
      //
      // ||grad^* u - grad u||_L^2(q)
      MPI_Reduce (&estgrad, &estGrad[adapt], 1, MPI_DOUBLE, MPI_SUM, 0, 
                  mpicomm);
      estGrad[adapt] = std::sqrt(estGrad[adapt]);
      
      // Number of mesh nodes at current step
      nnodes[adapt] = tmsh.num_global_nodes ();

      MPI_Barrier (mpicomm); if (rank == 0) { toc ("compute h and metrics "); }
      
      // Refinement procedure
      MPI_Barrier (mpicomm); if (rank == 0) { tic (); }

      // Break if the number of global nodes is too large
      if (tmsh.num_global_nodes () >= 10e7)
        break;
      else if (adapt < (adapt_refine_steps - 1))
        {   
          // Set marker for refinement
          tmsh.set_metrics_marker (estimator, tol*std::pow (.9, adapt),3,1,1);
          
          // Refine.
          tmsh.metrics_refine (1e5);
          std::cout << "tmsh.num_global_nodes ()= "
                    << tmsh.num_global_nodes ()
                    << std::endl;
          
          // Export new mesh
          tmsh.vtk_export ((std::string("p4est_adr_test_1_metrics_newmesh_")
                            + std::to_string(adapt)).c_str());

          // Break if the number of global nodes is too large
          if (tmsh.num_global_nodes () >= 10e7)
            {
              std::cout << "too many nodes!" << std::endl;
              break;
            }
        }
      MPI_Barrier (mpicomm); if (rank == 0) { toc ("refine "); }
    }
  
  if (rank == 0)
    for (unsigned step = 0; step < nnodes.size(); ++step)
      {
        std::cout << "\nStep " << step << ", #nodes: "
                  << nnodes[step] << ", h: "
                  << h_step[step] << std::endl;
        std::cout << "\n\tSolution estimator = " << estSol[step]
                  << "\n\tGradient estimator = " << estGrad[step]
                  << std::endl;
        std::cout << std::endl;                  
      }
  
  MPI_Finalize ();
  
  return 0;
}
