#include <bim_sparse_distributed.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
#include <cassert>
#include <limits>

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return 1; }



int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;

  using q1_vec = q1_vec<distributed_vector>;
  using gradient = gradient<distributed_vector>;

  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  // Create mesh
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  // Uniform refinement
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 2; ++cycle)
    {
      tmsh.refine (recursive, partforcoarsen);
      tmsh.set_refine_marker (uniform_refinement);
    }
  
  // Export initial mesh
  tmsh.vtk_export ("p4est_adr_test_2_distributed");
  
  // Mesh parameters
  std::vector<tmesh::idx_t> nnodes;     // number of nodes at every step
  std::vector<double> h_step;           // mesh size at every step
  constexpr unsigned refine_steps = 5;

  // Problem parameters
  constexpr double epsilon = 1e-6;      // diffusion coefficient
  constexpr double theta = M_PI / 4;    // angle

  // Adaptive refinement loop
  for (int adapt = 0; adapt < refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
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
                  
                  psi[quadrant->gt(ii)] = (std::cos(theta) * x +
                                           std::sin(theta) * y) / epsilon;
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

      bim2a_advection_diffusion (tmsh, alpha, psi, A);
      
      // Assemble right-hand side.
      q1_vec rhs(tmsh.num_owned_nodes ());
      bim2a_solution_with_ghosts(tmsh,rhs); 
      
      bim2a_rhs (tmsh, f, g, rhs);
      
      // Set boundary conditions.
      func u0  = [] (double x, double y) { return 0; };
      func u1  = [] (double x, double y) { return 1; };
      func u10 = [] (double x, double y) { return (y > 0.2) ? 0 : 1; };
      
      dirichlet_bcs bcs;
      bcs.push_back (std::make_tuple(0, 0, u10));
      bcs.push_back (std::make_tuple(0, 1, u0 ));
      bcs.push_back (std::make_tuple(0, 2, u1 ));
      bcs.push_back (std::make_tuple(0, 3, u0 ));
      
      bim2a_dirichlet_bc (tmsh, bcs, A, rhs);
      A.assemble();

      // Solve problem.
      std::cout << "Solving linear system.";
      
      mumps mumps_solver;
      
      std::vector<double> vals;
      std::vector<int> irow, jcol;
      
      A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
      mumps_solver.set_distributed_lhs_data (vals);

      mumps_solver.set_rhs_distributed (rhs);
      
      // Solve.
      mumps_solver.analyze ();
      mumps_solver.factorize ();
      mumps_solver.solve ();
      mumps_solver.cleanup ();

      // Get solution of linear system
      q1_vec mumps_result = mumps_solver.get_distributed_solution();

      q1_vec result (tmsh.num_owned_nodes());
      bim2a_solution_with_ghosts (tmsh, result);

      for (auto ii = result.get_range_start (); 
                ii != result.get_range_end (); 
                ++ii)
        result[ii] = mumps_result[ii];
      result.assemble (replace_op);

      // Export solution.
      tmsh.octbin_export ((std::string("p4est_adr_test_2_distributed_u_")
                           + std::to_string(adapt)).c_str(), result);

      std::cout << " Done." << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient and estimator.";
      
      gradient du = bim2c_quadtree_pde_recovered_gradient(tmsh, result);
      q2_vec u_star = bim2c_quadtree_pde_recovered_solution(tmsh,result,du);
      
      tmsh.octbin_export ((std::string("p4est_adr_test_2_distributed_du_x_")
                           + std::to_string(adapt)).c_str(), du.first);
      tmsh.octbin_export ((std::string("p4est_adr_test_2_distributed_du_y_")
                           + std::to_string(adapt)).c_str(), du.second);
      
      auto estimator = [& u_star, & result] (tmesh::quadrant_iterator q)
        { return estimator_sol (q, u_star, result); };
      
      double tol = 1e-6;
      tmsh.set_metrics_marker (estimator, tol, 4);
      
      // Compute metrics and h.
      std::vector<double> metrics(tmsh.num_local_quadrants ());
      
      double hx = 0, hy = 0,
        h = std::numeric_limits<double>::max (),
        global_h = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          metrics[quadrant->get_forest_quad_idx ()] =
            estimator (quadrant) * std::sqrt (tmsh.num_global_quadrants ())
            / tol;
          
          hx = quadrant->p(0, 1) - quadrant->p(0, 0);
          hy = quadrant->p(1, 2) - quadrant->p(1, 0);
          
          h = std::min(h, std::sqrt(hx*hx + hy*hy));
        }
      
      tmsh.octbin_export_quadrant((std::string(
                                    "p4est_adr_test_2_distributed_hx_")
                                    + std::to_string(adapt)).c_str(), metrics);
        
      MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      h_step.push_back (global_h);
      
      std::cout << " Done." << std::endl;
      
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Refine.
      tmsh.metrics_refine (1e5);
      
      tmsh.vtk_export ((std::string("p4est_adr_test_2_distributed_newmesh_")
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
