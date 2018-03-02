#include <bim_sparse.h>
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

static constexpr unsigned refine_steps = 20;

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
  
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 1; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  tmsh.vtk_export ("p4est_dr_test_4_marker");
  
  std::vector<tmesh::idx_t> nnodes;
  std::vector<double> h_step;
  std::vector<double> error;
  
  double delta1 = 1.5;
  double delta2 = 0.5;
  
  for (int adapt = 0; adapt < refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
      // Compute coefficients.
      double eps1 = 5e-5;
      double eps2 = 1e-1;
      
      auto c = [eps1, eps2] (double x)
        {
          return
          ((0.71875 - (0.125 * x - 0.375) * x) * eps2) /
          (
            (x - 1.5) * std::sqrt(eps1) *
            std::cosh((0.5 * x + 0.25) / std::sqrt(eps1)) - 
            2 * eps2 * std::sinh((0.5 * x + 0.25) / std::sqrt(eps1))
          );
        };
      
      auto d = [eps1, eps2] (double x)
        {
          return
          (
            (-1.8125 - (0.25 * x - 0.25) * x) * std::sqrt(eps1) *
            std::cosh((0.5 * x + 0.25) / std::sqrt(eps1)) +
            (x - 0.5) * eps2 * std::sinh((0.5 * x + 0.25) / std::sqrt(eps1))
          ) /
          (
            (x - 1.5) * std::sqrt(eps1) *
            std::cosh((0.5 * x + 0.25) / std::sqrt(eps1)) - 
            2 * eps2 * std::sinh((0.5 * x + 0.25) / std::sqrt(eps1))
          );
        };
      
      std::vector<double> alpha(tmsh.num_local_nodes (), eps1);
      std::vector<double> psi(tmsh.num_local_nodes (), 0);
      
      std::vector<double> delta(tmsh.num_local_quadrants (), 1);
      std::vector<double> zeta(tmsh.num_local_nodes (), 1);
      
      std::vector<double> f(tmsh.num_local_quadrants (), 1);
      std::vector<double> g(tmsh.num_local_nodes (), 1);
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int ii = 0; ii < 4; ++ii)
            if (quadrant->p(1, ii) > 0.5 * quadrant->p(0, ii) + 0.25)
              alpha[quadrant->t(ii)] = eps2;
          
          if (quadrant->centroid(1) > 0.5 * quadrant->centroid(0) + 0.25)
            {
              delta[quadrant->get_forest_quad_idx()] = 0;
              f[quadrant->get_forest_quad_idx()] = eps2;
            }
        }
      
      // Assemble system matrix and right-hand side.
      sparse_matrix A, M;
      A.resize(tmsh.num_global_nodes());
      M.resize(tmsh.num_global_nodes());
      bim2a_advection_eafe_diffusion (tmsh, alpha, psi, A);
      bim2a_reaction (tmsh, delta, zeta, M);
      A += M;
      
      std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      bim2a_rhs (tmsh, f, g, rhs);
      
      // Set boundary conditions.
      func u_ex =
        [eps1, eps2, c, d] (double x, double y)
          {
            if (y <= 0.5 * x + 0.25)
              return (1 + 2 * c(x) * std::sinh(y / std::sqrt(eps1)));
            else
              return (-0.5 * (y - 1) * (y + 2 * d(x)));
          };
      
      dirichlet_bcs bcs;
      bcs.push_back (std::make_tuple(0, 2, u_ex));
      bcs.push_back (std::make_tuple(0, 3, u_ex));
      
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
      tmsh.octbin_export ((std::string("p4est_dr_test_4_marker_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::vector<double> uex(tmsh.num_global_nodes(), 0);
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        for (int i = 0; i < 4; ++i)
          uex[quadrant->gt(i)] = u_ex(quadrant->p(0, i), quadrant->p(1, i));
      
      tmsh.octbin_export ((std::string("p4est_dr_test_4_marker_uex_")
                           + std::to_string(adapt)).c_str(), uex);
      
      std::cout << " Done." << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient, solution and estimator.";
      
      active_fun region0 = [] (tmesh::quadrant_iterator q)
        { return (q->centroid(1) <= 0.5 * q->centroid(0) + 0.25); };
      
      active_fun region1 = [] (tmesh::quadrant_iterator q)
        { return (q->centroid(1) > 0.5 * q->centroid(0) + 0.25); };
      
      gradient du0 = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs, region0);
      gradient du1 = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs, region1);
      
      q2_vec u_star0 = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du0);
      q2_vec u_star1 = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du1);
      
      tmsh.octbin_export ((std::string("p4est_dr_test_4_marker_du0_x_")
                           + std::to_string(adapt)).c_str(), du0.first);
      tmsh.octbin_export ((std::string("p4est_dr_test_4_marker_du0_y_")
                           + std::to_string(adapt)).c_str(), du0.second);
      
      tmsh.octbin_export ((std::string("p4est_dr_test_4_marker_du1_x_")
                           + std::to_string(adapt)).c_str(), du1.first);
      tmsh.octbin_export ((std::string("p4est_dr_test_4_marker_du1_y_")
                           + std::to_string(adapt)).c_str(), du1.second);
      
      auto refine_fun = [& delta1, & u_star0, & u_star1, & global_rhs, &tmsh] (tmesh::quadrant_iterator q)
        {
          if (q->centroid(1) <= 0.5 * q->centroid(0) + 0.25)
            return zz_marker_sol (q, u_star0, global_rhs,
                                  delta1 * 1e-10 / std::sqrt(tmsh.num_global_nodes()));
          else
            return zz_marker_sol (q, u_star1, global_rhs,
                                  delta1 * 1e-10 / std::sqrt(tmsh.num_global_nodes()));
        };
      
      auto coarsen_fun = [& delta2, & u_star0, & u_star1, & global_rhs, &tmsh] (tmesh::quadrant_iterator q)
        {
          if (q->centroid(1) <= 0.5 * q->centroid(0) + 0.25)
            return !zz_marker_sol (q, u_star0, global_rhs,
                                   delta2 * 1e-10 / std::sqrt(tmsh.num_global_nodes()));
          else
            return !zz_marker_sol (q, u_star1, global_rhs,
                                   delta2 * 1e-10 / std::sqrt(tmsh.num_global_nodes()));
        };
      
      // Compute h and error.
      double hx = 0, hy = 0,
             h = std::numeric_limits<double>::max (),
             global_h = 0;
      double err = 0, global_err = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          hx = quadrant->p(0, 1) - quadrant->p(0, 0);
          hy = quadrant->p(1, 2) - quadrant->p(1, 0);
          
          h = std::min(h, std::sqrt(hx*hx + hy*hy));
          
          err += std::pow(l2_error(quadrant, u_ex, global_rhs), 2);
        }
      
      MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      MPI_Reduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      global_err = std::sqrt(global_err);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      h_step.push_back (global_h);
      error.push_back (global_err);
      
      std::cout << " Done." << std::endl;
      
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Coarsen and refine.
      tmsh.set_coarsen_marker (coarsen_fun);
      tmsh.set_refine_marker (refine_fun);
      
      tmsh.coarsen (recursive, partforcoarsen, 0);
      tmsh.refine (recursive, partforcoarsen);
      
      tmsh.vtk_export ((std::string("p4est_dr_test_4_marker_newmesh_")
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
