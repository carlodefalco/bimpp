#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <vector>
#include <cassert>

// Define mesh.
constexpr p4est_topidx_t simple_conn_num_vertices = 6;
constexpr p4est_topidx_t simple_conn_num_trees = 2;
const double simple_conn_p[simple_conn_num_vertices*2] = {0, 0, 1, 0, 1, 0.5, 0, 0.5, 1, 1, 0, 1};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = {1, 2, 3, 4, 1, 4, 3, 5, 6, 1};

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return 1; }

static constexpr unsigned refine_steps = 10;

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  tmsh.set_replace_fun (tmesh::userint_replace);
  
  tmsh.set_refine_marker (uniform_refinement);
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 2; ++cycle)
    tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_dr_test_3_metrics");
  
  std::vector<tmesh::idx_t> nnodes;
  
  for (int adapt = 0; adapt < refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
      // Compute coefficients.
      double eps1 = 5e-7;
      double eps2 = 1e-6;
      
      double c = -0.4375 * eps2 /
        (0.5 * std::sqrt(eps1) * std::cosh(0.5 / std::sqrt(eps1)) +
         eps2 * std::sinh(0.5 / std::sqrt(eps1)));
      
      double d = 1.75 * std::sqrt(eps1) * std::cosh(0.5 / std::sqrt(eps1)) /
        (std::sqrt(eps1) * std::cosh(0.5 / std::sqrt(eps1)) +
         2 * eps2 * std::sinh(0.5 / std::sqrt(eps1)));
      
      std::vector<double> alpha(tmsh.num_local_quadrants (), eps1);
      std::vector<double> psi(tmsh.num_local_nodes (), 0);
      
      std::vector<double> delta(tmsh.num_local_quadrants (), 1);
      std::vector<double> zeta(tmsh.num_local_nodes (), 1);
      
      std::vector<double> f(tmsh.num_local_quadrants (), 1);
      std::vector<double> g(tmsh.num_local_nodes (), 1);
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          if (quadrant->p(1, 0) >= 0.5)
            {
              alpha[quadrant->get_forest_quad_idx()] = eps2;
              delta[quadrant->get_forest_quad_idx()] = 0;
              f[quadrant->get_forest_quad_idx()] = eps2;
            }
        }
      
      // Assemble system matrix and right-hand side.
      sparse_matrix A, M;
      A.resize(tmsh.num_global_nodes());
      M.resize(tmsh.num_global_nodes());
      bim2a_advection_diffusion (tmsh, alpha, psi, A);
      bim2a_reaction (tmsh, delta, zeta, M);
      A += M;
      
      std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      bim2a_rhs (tmsh, f, g, rhs);
      
      // Set boundary conditions.
      func u_ex =
        [eps1, eps2, c, d] (double x, double y)
          {
            if (y <= 0.5)
              return (1 + 2 * c * std::sinh(y / std::sqrt(eps1)));
            else
              return (-0.5 * (y - 1) * (y + 2 * d));
          };
      
      dirichlet_bcs bcs;
      bcs.push_back (std::make_tuple(0, 0, u_ex));
      bcs.push_back (std::make_tuple(1, 1, u_ex));
      
      bim2a_dirichlet_bc (tmsh, bcs, A, rhs);
      
      // Solve problem.
      std::cout << "Solving linear system.";
      
      mumps mumps_solver;
      
      std::vector<double> vals;
      std::vector<int> irow, jcol;
      
      A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
      mumps_solver.set_distributed_lhs_data (vals);
      
      // Reduce rhs (so that rank 0 has the actual rhs).
      std::vector<double> global_rhs(tmsh.num_global_nodes(), 0);
      MPI_Reduce(rhs.data(), global_rhs.data(), rhs.size(),
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      
      if (rank == 0)
        mumps_solver.set_rhs (global_rhs);
      
      // Solve.
      mumps_solver.analyze ();
      mumps_solver.factorize ();
      mumps_solver.solve ();
      mumps_solver.cleanup ();
      
      // Export solution.
      MPI_Bcast(global_rhs.data(), global_rhs.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      tmsh.octbin_export ((std::string("p4est_dr_test_3_metrics_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::vector<double> uex(tmsh.num_global_nodes(), 0);
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        for (int i = 0; i < 4; ++i)
          uex[quadrant->gt(i)] = u_ex(quadrant->p(0, i), quadrant->p(1, i));
      
      tmsh.octbin_export ((std::string("p4est_dr_test_3_metrics_uex_")
                           + std::to_string(adapt)).c_str(), uex);
      
      std::cout << " Done." << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient and estimator.";
      
      gradient du = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs);
      q2_vec u_star = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du);
      
      tmsh.octbin_export ((std::string("p4est_dr_test_3_metrics_du_x_")
                           + std::to_string(adapt)).c_str(), du.first);
      tmsh.octbin_export ((std::string("p4est_dr_test_3_metrics_du_y_")
                           + std::to_string(adapt)).c_str(), du.second);
      
      auto estimator = [& u_star, & global_rhs] (tmesh::quadrant_iterator q)
        { return estimator_sol (q, u_star, global_rhs); };
      
      nnodes.push_back (tmsh.num_global_nodes ());
      
      std::cout << " Done." << std::endl;
      
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Refine.
      tmsh.set_metrics_marker (estimator, 1e-10, 4);
      tmsh.metrics_refine ();
      
      tmsh.vtk_export ((std::string("p4est_dr_test_3_metrics_newmesh_")
                        + std::to_string(adapt)).c_str());
    }
  
  if (rank == 0)
    for (unsigned step = 0; step < nnodes.size(); ++step)
      std::cout << "Step " << step << ", #nodes: "
                << nnodes[step] << std::endl;
  
  MPI_Finalize ();
  
  return 0;
}
