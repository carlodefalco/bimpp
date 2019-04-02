#include <mumps_class.h>
#include <quad_operators_3d.h>
#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
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
  
  using q1_vec          = q1_vec<distributed_vector>;
  using gradient3       = gradient3<std::vector<double>>;
  using idx_t           = tmesh_3d::idx_t;

  // Compute radius in the given point
  func3 rho2 = [] (double x, double y, double z) -> double
    {
      return ((x - 0.5) * (x - 0.5) + 
              (y - 0.5) * (y - 0.5) +
              (z - 0.5) * (z - 0.5));
    };
    
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  // Number of refinement steps 
  constexpr unsigned unif_refine_steps = 2;   // initial uniform refinement 
  constexpr unsigned adapt_refine_steps = 5;  // adaptive refinement 

  // Mesh parameters
  std::vector<idx_t> nnodes;            // number of nodes at every step
  std::vector<double> h_step;           // mesh size at every step

  // Error at every step
  std::vector<double> error;

  // Problem parameters
  constexpr double inv_epsilon = 1e11;  // 1 / epsilon
  constexpr double R = 0.25;            // radius of the internal sphere
  constexpr double kS = 1;            	// diffusion coefficient in the sphere
  constexpr double kG = 1;              // diffusion coefficient outside

  // Exact solution
  func3 u_ex = [rho2,R] (double x, double y, double z) 
  {
    double R2 = R*R;
    double r2 = rho2 (x, y, z);
    return (r2 > R2 ?  std::sin (r2) : std::sin (R2));
  }; 

  // Load term
  func3 load = [rho2,R,inv_epsilon] (double x, double y, double z)
  {
    double R2 = R*R;
    double r2 = rho2 (x, y, z);
    return (r2 > R2 ?
            (6 * std::cos (r2) - 4 * r2 * std::sin (r2)) :
            inv_epsilon * std::sin (R2));
  };

  // Diffusion term
  func3 diffusion = [rho2,R,kG,kS] (double x, double y, double z)
  {
    double R2 = R*R;
    double r2 = rho2 (x, y, z);
    return (r2 > R2 ? kG : kS);
  }; 

  // Reaction term
  func3 reaction = [rho2,R,inv_epsilon] (double x, double y, double z)
  {
    double R2 = R*R;
    double r2 = rho2 (x, y, z);
    return (r2 > R2 ? 0.0 : inv_epsilon);
  };
  
  // Mesh generation
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  // Initial uniform refinement
  recursive = 0; partforcoarsen = 1;
  for (unsigned cycle = 0; cycle < unif_refine_steps; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }

  // Export initial mesh
  tmsh.vtk_export ("p4est_dr_test_2_metrics_initial_mesh");
  
  // Adaptive refinement loop
  for (unsigned adapt = 0; adapt < adapt_refine_steps; ++adapt)
    {
      std::cout << "*** Step "
                << adapt << " (rank "
                << rank << ") ***" << std::endl;
      
      // Compute coefficients   
      //
      // diffusion   
      std::vector<double> alpha (tmsh.num_local_quadrants (), 0.);
      q1_vec psi (tmsh.num_owned_nodes ());
      //
      // reaction
      std::vector<double> delta (tmsh.num_local_quadrants (), -1.);
      q1_vec zeta (tmsh.num_owned_nodes ());
      //
      // rhs
      std::vector<double> f (tmsh.num_local_quadrants (), -1.);
      q1_vec g (tmsh.num_owned_nodes ());
      //
      double x = .0, y = .0, z = .0;
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int ii = 0; ii < 8; ++ii)
            {
              x = quadrant->p(0,ii);
              y = quadrant->p(1,ii);
              z = quadrant->p(2,ii);

              if (! quadrant->is_hanging(ii))
                {
                  psi[quadrant->gt(ii)] = 0.;
                  zeta[quadrant->gt(ii)] = reaction (x, y, z);
                  g[quadrant->gt(ii)] = load (x, y, z);
                }
              else
                {
                  psi[quadrant->gt(ii)] += 0.;
                  zeta[quadrant->gt(ii)] += 0.;
                  g[quadrant->gt(ii)] += 0.;
                }
            }
          alpha[quadrant->get_forest_quad_idx()] = diffusion(x, y, z);
        }
      psi.assemble (replace_op);
      zeta.assemble (replace_op);
      g.assemble (replace_op);
      
      // Assemble system matrix and right-hand side.
      distributed_sparse_matrix A;
      A.set_ranges(tmsh.num_owned_nodes());
      //
      // advection_diffusion
      bim3a_advection_diffusion (tmsh, alpha, psi, A);
      //
      // reaction
      bim3a_reaction (tmsh, delta, zeta, A);
      
      // rhs
      q1_vec rhs (tmsh.num_owned_nodes ());
      bim3a_solution_with_ghosts (tmsh, rhs);
      //
      bim3a_rhs (tmsh, f, g, rhs);

      /// CDF : end checked

      // Set boundary conditions.      
      dirichlet_bcs3 bcs;
      for (int i = 0; i < 6; ++i)
        bcs.push_back (std::make_tuple(0, i, u_ex));
      
      bim3a_dirichlet_bc (tmsh, bcs, A, rhs);
      
      // Solve problem.
      MPI_Barrier(mpicomm);
      std::cout << "Solving linear system. (rank " << rank << ")" << std::endl;
      
      // Initialize MUMPS solver
      mumps mumps_solver;//(true);
      
      // Set distributed structure of lhs
      std::vector<double> vals;
      std::vector<int> irow, jcol;
      
      A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
      mumps_solver.set_distributed_lhs_data (vals);
      
      // Set distributed structure of rhs
      mumps_solver.set_rhs_distributed (rhs);

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
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      // Compute exact solution on the mesh
      std::vector<double> uex(tmsh.num_global_nodes(), 0);
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        for (int i = 0; i < 8; ++i)
          if (! quadrant->is_hanging (i))
            uex[quadrant->gt(i)] = u_ex(quadrant->p(0, i), 
                                        quadrant->p(1, i),
                                        quadrant->p(2, i));
      
      // Export exact solution
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_uex_")
                           + std::to_string(adapt)).c_str(), uex);
      
      std::cout << "Done. (rank " << rank << ")" << std::endl;
      
      // Compute reconstructed gradient and solution
      std::cout << "Computing reconstructed gradient and estimator.";
      
      // Activation function outside the sphere
      active_fun3 regionG = [rho2,R] (tmesh_3d::quadrant_iterator q)
      { 
        return (rho2(q->centroid(0),q->centroid(1),q->centroid(2)) > (R*R));
      };
      
      // Activation function inside the sphere
      active_fun3 regionS = [rho2,R] (tmesh_3d::quadrant_iterator q)
      { 
        return (rho2(q->centroid(0),q->centroid(1),q->centroid(2)) <= (R*R));
      };
      
      // Gradient outside the sphere
      std::cout << "\tgradient (rank " << rank << ")" << std::endl;
      gradient3 du0 = bim3c_quadtree_pde_recovered_gradient(tmsh, 
                                                            global_rhs, 
                                                            regionG);
      // Gradient inside the sphere
      gradient3 du1 = bim3c_quadtree_pde_recovered_gradient(tmsh, 
                                                            global_rhs, 
                                                            regionS);
      // Reconstructed solution outside the sphere
      std::cout << "\tsolution (rank " << rank << ")" << std::endl;
      q2_vec3 u_star0 = bim3c_quadtree_pde_recovered_solution(tmsh, 
                                                              global_rhs,
                                                              du0);
      // Reconstructed solution inside the sphere
      q2_vec3 u_star1 = bim3c_quadtree_pde_recovered_solution(tmsh, 
                                                              global_rhs, 
                                                              du1);
      
      // Export reconstructed gradients:
      //
      // derivative along x outside the sphere 
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_du0_x_")
                           + std::to_string(adapt)).c_str(), std::get<0>(du0));
      //
      // derivative along y outside the sphere 
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_du0_y_")
                           + std::to_string(adapt)).c_str(), std::get<1>(du0));
      //
      // derivative along z outside the sphere 
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_du0_z_")
                           + std::to_string(adapt)).c_str(), std::get<2>(du0));
      //
      // derivative along x inside the sphere 
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_du1_x_")
                           + std::to_string(adapt)).c_str(), std::get<0>(du1));
      //
      // derivative along y inside the sphere 
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_du1_y_")
                           + std::to_string(adapt)).c_str(), std::get<1>(du1));
      //
      // derivative along z inside the sphere 
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_du1_z_")
                           + std::to_string(adapt)).c_str(), std::get<2>(du1));
      
      // Solution estimator
      auto estimator = [rho2, &u_star0, &u_star1, &global_rhs, R] 
                                                (tmesh_3d::quadrant_iterator q)
        {
          if (rho2(q->centroid(0),q->centroid(1),q->centroid(2)) > (R*R))
            return estimator_sol (q, u_star0, global_rhs);
          else
            return estimator_sol (q, u_star1, global_rhs);
        };
      
      std::cout << " Done. (rank " << rank << ")" << std::endl;
      
      // Compute h and error.
      std::cout << "\tmetrics (rank " << rank << ")" << std::endl;
      double  hx = 0, hy = 0, hz = 0,
              h = std::numeric_limits<double>::max (),
              global_h = 0;
      double err = 0, global_err = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          hx = quadrant->p(0, 1) - quadrant->p(0, 0);
          hy = quadrant->p(1, 7) - quadrant->p(1, 0);
          hy = quadrant->p(2, 7) - quadrant->p(2, 0);
          
          h = std::min(h, std::sqrt(hx*hx + hy*hy + hz*hz));
          
          // ||u - u_ex||_L^2(q)
          err += std::pow(l2_error(quadrant, u_ex, global_rhs), 2);
        }
      
      // Global mesh size and global error
      MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      MPI_Reduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      global_err = std::sqrt(global_err);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      h_step.push_back (global_h);
      error.push_back (global_err);

      std::cout << " Done. (rank " << rank << ")\n" << std::endl;
      
      // Break if the number of global nodes is too large
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Refine.
      tmsh.set_metrics_marker (estimator, 1e-3, 4);
      tmsh.metrics_refine (1e3);
      
      // Export new mesh
      tmsh.vtk_export ((std::string("p4est_dr_test_2_metrics_newmesh_")
                        + std::to_string(adapt)).c_str());
    }
  
  if (rank == 0)
    for (unsigned step = 0; step < nnodes.size(); ++step)
      std::cout << "Step " << step << ", #nodes: "
                << nnodes[step] << ", h: "
                << h_step[step] << ", error: "
                << error[step] << std::endl;
  
  MPI_Finalize ();
  
  return 0;
}