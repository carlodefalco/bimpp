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
  
  tmsh.set_replace_fun (tmesh::user_int_replace);
  
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 2; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  tmsh.vtk_export ("p4est_dr_test_2_metrics");
  
  std::vector<tmesh::idx_t> nnodes;
  std::vector<double> h_step;
  
  for (int adapt = 0; adapt < refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
      // Assemble matrix and right-hand side.
      sparse_matrix A, M;
      A.resize(tmsh.num_global_nodes());
      M.resize(tmsh.num_global_nodes());
      
      double epsilon = std::pow(2, -10);
      std::vector<double> alpha(tmsh.num_local_quadrants (), epsilon);
      std::vector<double> psi(tmsh.num_global_nodes (), 0);
      
      std::vector<double> delta(tmsh.num_local_quadrants (), 1);
      std::vector<double> zeta(tmsh.num_global_nodes (), 1);
      
      std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      
      std::vector<double> f(tmsh.num_local_quadrants (), 1);
      std::vector<double> g(tmsh.num_global_nodes (), 0);
      
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
                  
                  zeta[quadrant->gt(ii)] = 1 + x * x * y * y;
                  g   [quadrant->gt(ii)] = 1 + 2 * x * y;
                }
            }
        }
        
      // Reduce coefficients.
      std::vector<double> global_zeta(tmsh.num_global_nodes(), 0);
      MPI_Allreduce(zeta.data(), global_zeta.data(), zeta.size(),
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
      std::vector<double> global_g(tmsh.num_global_nodes(), 0);
      MPI_Allreduce(g.data(), global_g.data(), g.size(),
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
      bim2a_advection_diffusion (tmsh, alpha, psi, A);
      bim2a_reaction (tmsh, delta, global_zeta, M);
      A += M;
      
      bim2a_rhs (tmsh, f, global_g, rhs);
      
      // Set boundary conditions.
      func g1 = [] (double x, double y) { return 1; };
      func g2 = [] (double x, double y) { return 1 - y * y; };
      func g3 = [] (double x, double y) { return 1 - x * x; };
      
      dirichlet_bcs bcs;
      bcs.push_back (std::make_tuple(0, 0, g1));
      bcs.push_back (std::make_tuple(0, 1, g2));
      bcs.push_back (std::make_tuple(0, 2, g1));
      bcs.push_back (std::make_tuple(0, 3, g3));
      
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
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::cout << " Done." << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient and estimator.";
      
      gradient du = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs);
      q2_vec u_star = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du);
      
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_du_x_")
                           + std::to_string(adapt)).c_str(), du.first);
      tmsh.octbin_export ((std::string("p4est_dr_test_2_metrics_du_y_")
                           + std::to_string(adapt)).c_str(), du.second);
      
      auto estimator = [& du, & global_rhs] (tmesh::quadrant_iterator q)
        { return estimator_grad (q, du, global_rhs); };
      
      double tol = 1e-3;
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
      
      tmsh.octbin_export_quadrant ((std::string("p4est_dr_test_2_metrics_hx_")
                                   + std::to_string(adapt)).c_str(), metrics);
      
      MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      h_step.push_back (global_h);
      
      std::cout << " Done." << std::endl;
      
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Refine.
      tmsh.metrics_refine (1e5);
      
      tmsh.vtk_export ((std::string("p4est_dr_test_2_metrics_newmesh_")
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
