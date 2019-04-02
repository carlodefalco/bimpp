#include <bim_sparse_distributed.h>
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
  using gradient3 = gradient3<std::vector<double>>;
  using idx_t = tmesh_3d::idx_t;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  // Number of refinement steps
  constexpr unsigned unif_refine_steps = 3;   // initial uniform refinement
  constexpr unsigned adapt_refine_steps = 4;  // adaptive refinement
  
  // Mesh parameters
  std::vector<idx_t>    nnodes;         // number of nodes at every step
  std::vector<double>   h_step;         // mesh size at every step

  // Error at every step
  std::vector<double> error;

  // Problem parameters
  constexpr double epsilon = 1e-6;              // diffusion coefficient
  double n_coeff = 1/std::sqrt(3.0);            // normalization coefficient

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
      std::cout << "*** Step " 
                << adapt << " (rank " 
                << rank << ") ***" << std::endl;
      
      // Compute coefficients
      //
      // diffusion
      std::vector<double> alpha(tmsh.num_local_quadrants (), epsilon);
      q1_vec psi(tmsh.num_owned_nodes ());
      //
      // rhs
      std::vector<double> f(tmsh.num_local_quadrants (), 0);
      q1_vec g(tmsh.num_owned_nodes ());
      //
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
                  
                  // CCI: divisione per epsilon --> matrice singolare (?)
                  psi[quadrant->gt(ii)] = (x + y - z) * n_coeff;// / epsilon;

                  g[quadrant->gt(ii)] = 0.;
                }
              else
                {
                  psi[quadrant->gt(ii)] += 0.;
                  g[quadrant->gt(ii)] += 0.;
                }
            }
        }
      psi.assemble(max_op);
      g.assemble(replace_op);

      // Assemble matrix.
      distributed_sparse_matrix A;
      A.set_ranges(tmsh.num_owned_nodes());
      //
      // advection_diffusion
      bim3a_advection_diffusion (tmsh, alpha, psi, A);
      
      // Assemble right-hand side.
      q1_vec rhs(tmsh.num_owned_nodes ());
      bim3a_solution_with_ghosts(tmsh,rhs);  
      
      bim3a_rhs (tmsh, f, g, rhs);

      // Set boundary conditions.
      func3 u0  = [] (double x, double y, double z) { return 0; };
      func3 u10 = [] (double x, double y, double z) 
      {
        double d = 1.;
        if (x <= 0.5 && y <= (- x + 0.5))
          d = 0.;
        return d; 
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
      
      // Solve problem.
      std::cout << "Solving linear system (rank " << rank << ")" << std::endl;

      // Initialize MUMPS solver
      mumps mumps_solver(true);
      
      // Set distributed structure of lhs
      std::vector<double> vals;
      std::vector<int> irow, jcol;
      
      A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
      mumps_solver.set_distributed_lhs_data (vals);
      
      // Set distributed structure of rhs
      mumps_solver.set_rhs_distributed(rhs);

      // Solve.
      std::cout << "\tanalyze (rank " << rank << ")" << std::endl;
      mumps_solver.analyze ();
      std::cout << "\tfactorize (rank " << rank << ")" << std::endl;
			mumps_solver.factorize ();
      std::cout << "\tsolve (rank " << rank << ")" << std::endl;
			mumps_solver.solve ();
      std::cout << "\tcleanup (rank " << rank << ")" << std::endl;
      mumps_solver.cleanup ();

      // Get solution on rank 0...
      q1_vec rhs_on_0 = mumps_solver.get_distributed_solution();
      //
      // ...and send it to all ranks
      unsigned size_global_rhs = 0;
      std::vector<double> global_rhs;
      if (rank == 0)
        {
          global_rhs = rhs_on_0.get_owned_data();
          size_global_rhs = global_rhs.size();
        }
      //
      MPI_Bcast(&size_global_rhs, 1, MPI_UNSIGNED, 0 , mpicomm);
      //
      if (rank != 0)
        global_rhs.resize(size_global_rhs);
      //
      MPI_Bcast(global_rhs.data(), size_global_rhs, MPI_DOUBLE, 0, mpicomm);

      // Export solution.
      tmsh.octbin_export ((std::string("p4est_adr_test_1_metrics_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::cout << "Done (rank " << rank << ")" << std::endl;
      
      // Compute reconstructed gradient and reconstructed solution
      std::cout << "Computing reconstructed gradient and estimator (rank "
      					<< rank << ")" << std::endl;
      
      std::cout << "\tgradient (rank " << rank << ")" << std::endl;
      gradient3 du = bim3c_quadtree_pde_recovered_gradient(tmsh, global_rhs);

      std::cout << "\tsolution (rank " << rank << ")" << std::endl;
      q2_vec3 u_star = bim3c_quadtree_pde_recovered_solution(tmsh,
                                                              global_rhs,
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
