#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
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
  
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 2; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  tmsh.vtk_export ("p4est_dr_test_5_uniform");
  
  std::vector<tmesh::idx_t> nnodes;
  std::vector<double> h_step;
  std::vector<double> error, error_du;
  
  for (int adapt = 0; adapt < refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
      // Compute coefficients.
      double R = 0.25;
      
      double kG = 1;
      double kS = 100;
      
      std::vector<double> alpha(tmsh.num_global_nodes (), kG);
      std::vector<double> psi(tmsh.num_global_nodes (), 0);
      
      std::vector<double> f(tmsh.num_local_quadrants (), 1);
      std::vector<double> g(tmsh.num_global_nodes (), 1);
      
      double x = 0, y = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int ii = 0; ii < 4; ++ii)
            {
              x = quadrant->p(0, ii);
              y = quadrant->p(1, ii);
                
              if (! quadrant->is_hanging (ii) &&
                  std::pow(x - 0.5, 2) +
                  std::pow(y - 0.5, 2) <=
                  std::pow(R, 2))
                alpha[quadrant->gt(ii)] = kS;
            }
        }
      
      // Assemble system matrix and right-hand side.
      sparse_matrix A;
      A.resize(tmsh.num_global_nodes());
      
      // Reduce coefficients.
      std::vector<double> global_alpha(tmsh.num_global_nodes(), 0);
      MPI_Allreduce(alpha.data(), global_alpha.data(), alpha.size(),
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
      bim2a_advection_eafe_diffusion (tmsh, global_alpha, psi, A);
      
      std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      bim2a_rhs (tmsh, f, g, rhs);
      
      // Set boundary conditions.
      func u_ex =
        [kG, kS, R] (double x, double y)
        {
          if (std::pow(x - 0.5, 2) +
              std::pow(y - 0.5, 2) >
              std::pow(R, 2))
            return
              (1.0 / 8 - 1 / (4 * kG) *
               (std::pow(x - 0.5, 2) + std::pow(y - 0.5, 2)));
          else
            return
              (1.0 / 8 - 1 / (4 * kS) *
               (std::pow(x - 0.5, 2) + std::pow(y - 0.5, 2)) -
               std::pow(R, 2) / 4 * (1 - 1 / kS));
        };
      
      dirichlet_bcs bcs;
      for (int i = 0; i < 4; ++i)
        bcs.push_back (std::make_tuple(0, i, u_ex));
      
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
      tmsh.octbin_export ((std::string("p4est_dr_test_5_uniform_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::vector<double> uex(tmsh.num_global_nodes(), 0);
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        for (int i = 0; i < 4; ++i)
          uex[quadrant->gt(i)] = u_ex(quadrant->p(0, i), quadrant->p(1, i));
      
      tmsh.octbin_export ((std::string("p4est_dr_test_5_uniform_uex_")
                           + std::to_string(adapt)).c_str(), uex);
      
      std::cout << " Done." << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient and solution.";
      
      active_fun regionG = [R] (tmesh::quadrant_iterator q)
        { return (std::pow(q->centroid(0) - 0.5, 2) +
                  std::pow(q->centroid(1) - 0.5, 2) >
                  std::pow(R, 2)); };
      
      active_fun regionS = [R] (tmesh::quadrant_iterator q)
        { return (std::pow(q->centroid(0) - 0.5, 2) +
                  std::pow(q->centroid(1) - 0.5, 2) <=
                  std::pow(R, 2)); };
      
      gradient<std::vector<double>> du0 = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs, regionG);
      gradient<std::vector<double>> du1 = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs, regionS);
      
      q2_vec u_star0 = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du0);
      q2_vec u_star1 = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du1);
      
      tmsh.octbin_export ((std::string("p4est_dr_test_5_uniform_du0_x_")
                           + std::to_string(adapt)).c_str(), du0.first);
      tmsh.octbin_export ((std::string("p4est_dr_test_5_uniform_du0_y_")
                           + std::to_string(adapt)).c_str(), du0.second);
      
      tmsh.octbin_export ((std::string("p4est_dr_test_5_uniform_du1_x_")
                           + std::to_string(adapt)).c_str(), du1.first);
      tmsh.octbin_export ((std::string("p4est_dr_test_5_uniform_du1_y_")
                           + std::to_string(adapt)).c_str(), du1.second);
      
      // Compute h and error.
      double hx = 0, hy = 0,
        h = std::numeric_limits<double>::max (),
        global_h = 0;
      
      func du_x_ex =
        [kG, kS, R] (double x, double y)
        {
          if (std::pow(x - 0.5, 2) +
              std::pow(y - 0.5, 2) >
              std::pow(R, 2))
            return (-(x - 0.5) / (2 * kG));
          else
            return (-(x - 0.5) / (2 * kS));
        };
      
      func du_y_ex =
        [kG, kS, R] (double x, double y)
        {
          if (std::pow(x - 0.5, 2) +
              std::pow(y - 0.5, 2) >
              std::pow(R, 2))
            return (-(y - 0.5) / (2 * kG));
          else
            return (-(y - 0.5) / (2 * kS));
        };
          
      double err = 0, global_err = 0;
      double err_du = 0, global_err_du = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          hx = quadrant->p(0, 1) - quadrant->p(0, 0);
          hy = quadrant->p(1, 2) - quadrant->p(1, 0);
          
          h = std::min(h, std::sqrt(hx*hx + hy*hy));
          
          err += std::pow(l2_error(quadrant, u_ex, global_rhs), 2);
          
          if (regionG(quadrant))
            err_du += std::pow(semih1_star_error(quadrant, du_x_ex, du_y_ex, du0), 2);
          else
            err_du += std::pow(semih1_star_error(quadrant, du_x_ex, du_y_ex, du0), 2);
        }
      
      MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      MPI_Reduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      MPI_Reduce(&err_du, &global_err_du, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      
      global_err = std::sqrt(global_err);
      global_err_du = std::sqrt(global_err_du);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      h_step.push_back (global_h);
      error.push_back (global_err);
      error_du.push_back (global_err_du);
      
      std::cout << " Done." << std::endl;
      
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Refine.
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
      
      tmsh.vtk_export ((std::string("p4est_dr_test_5_uniform_newmesh_")
                        + std::to_string(adapt)).c_str());
    }
  
  if (rank == 0)
    for (unsigned step = 0; step < nnodes.size(); ++step)
      std::cout << "Step " << step << ", #nodes: "
                << nnodes[step] << ", h: "
                << h_step[step] << ", error: "
                << error[step] <<  ", error_du^*: "
                << error_du[step] << std::endl;
  
  MPI_Finalize ();
  
  return 0;
}
