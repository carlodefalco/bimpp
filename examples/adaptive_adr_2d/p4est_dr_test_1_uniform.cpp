#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
#include <cassert>

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
  for (int cycle = 0; cycle < 2; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  tmsh.vtk_export ("p4est_dr_test_1_uniform");
  
  std::vector<tmesh::idx_t> nnodes;
  std::vector<double> error, error_du_x, error_du_y;
  
  double delta1 = 1.5;
  double delta2 = 0.5;
  
  for (int adapt = 0; adapt < refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
      // Assemble matrix.
      sparse_matrix A, M;
      A.resize(tmsh.num_global_nodes());
      M.resize(tmsh.num_global_nodes());
      
      double epsilon = 1e-5;
      std::vector<double> alpha(tmsh.num_local_quadrants (), epsilon);
      std::vector<double> psi(tmsh.num_local_nodes (), 0);
      
      std::vector<double> delta(tmsh.num_local_quadrants (), 1);
      std::vector<double> zeta(tmsh.num_local_nodes (), 1);
      
      bim2a_advection_diffusion (tmsh, alpha, psi, A);
      bim2a_reaction (tmsh, delta, zeta, M);
      A += M;
      
      // Assemble right-hand side.
      std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      
      std::vector<double> f(tmsh.num_local_quadrants (), 1);
      std::vector<double> g(tmsh.num_local_nodes (), 0);
      
      double x = 0, y = 0;
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int ii = 0; ii < 4; ++ii)
            {
              x = quadrant->p(0, ii);
              y = quadrant->p(1, ii);
              
              g[quadrant->t(ii)] = 1 - std::sinh(x / std::sqrt(epsilon)) *
                                       std::sinh(y / std::sqrt(epsilon)) /
                                       std::sinh(1 / std::sqrt(epsilon)) /
                                       std::sinh(1 / std::sqrt(epsilon));
            }
        }
        
      bim2a_rhs (tmsh, f, g, rhs);
      
      // Set boundary conditions.
      func u_ex =
        [epsilon] (double x, double y)
        { return (1 - std::sinh(x / std::sqrt(epsilon)) / std::sinh(1 / std::sqrt(epsilon))) *
                 (1 - std::sinh(y / std::sqrt(epsilon)) / std::sinh(1 / std::sqrt(epsilon))); };
                 
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
      tmsh.octbin_export ((std::string("p4est_dr_test_1_uniform_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::vector<double> uex(tmsh.num_global_nodes(), 0);
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        for (int i = 0; i < 4; ++i)
          uex[quadrant->gt(i)] = u_ex(quadrant->p(0, i), quadrant->p(1, i));
      
      tmsh.octbin_export ((std::string("p4est_dr_test_1_uniform_uex_")
                           + std::to_string(adapt)).c_str(), uex);
      
      std::cout << " Done." << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient and solution.";
      
      gradient du = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs);
      q2_vec u_star = bim2c_quadtree_pde_recovered_solution(tmsh, global_rhs, du);
      
      tmsh.octbin_export ((std::string("p4est_dr_test_1_uniform_du_x_")
                           + std::to_string(adapt)).c_str(), du.first);
      tmsh.octbin_export ((std::string("p4est_dr_test_1_uniform_du_y_")
                           + std::to_string(adapt)).c_str(), du.second);
      
      // Compute error.
      func du_x_ex =
        [epsilon] (double x, double y)
        { return (- std::cosh(x / std::sqrt(epsilon)) / std::sinh(1 / std::sqrt(epsilon))) *
                 (1 - std::sinh(y / std::sqrt(epsilon)) / std::sinh(1 / std::sqrt(epsilon))) / std::sqrt(epsilon); };
      
      func du_y_ex =
        [epsilon] (double x, double y)
        { return (1 - std::sinh(x / std::sqrt(epsilon)) / std::sinh(1 / std::sqrt(epsilon))) *
                 (- std::cosh(y / std::sqrt(epsilon)) / std::sinh(1 / std::sqrt(epsilon))) / std::sqrt(epsilon); };
      
      double err = 0, global_err = 0;
      double err_du_x = 0, global_err_du_x = 0;
      double err_du_y = 0, global_err_du_y = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          err += std::pow(l2_error(quadrant, u_ex, global_rhs), 2);
          err_du_x += std::pow(l2_error(quadrant, du_x_ex, du.first), 2);
          err_du_y += std::pow(l2_error(quadrant, du_y_ex, du.second), 2);
        }
      
      MPI_Reduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      MPI_Reduce(&err_du_x, &global_err_du_x, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      MPI_Reduce(&err_du_y, &global_err_du_y, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      
      global_err = std::sqrt(global_err);
      global_err_du_x = std::sqrt(global_err_du_x);
      global_err_du_y = std::sqrt(global_err_du_y);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      error.push_back (global_err);
      error_du_x.push_back (global_err_du_x);
      error_du_y.push_back (global_err_du_y);
      
      std::cout << " Done." << std::endl;
      
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Refine.
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
      
      tmsh.vtk_export ((std::string("p4est_dr_test_1_uniform_newmesh_")
                        + std::to_string(adapt)).c_str());
    }
  
  if (rank == 0)
    for (unsigned step = 0; step < nnodes.size(); ++step)
      std::cout << "Step " << step << ", #nodes: "
                << nnodes[step] << ", error: "
                << error[step] <<  ", error_du_x: "
                << error_du_x[step] <<  ", error_du_y: "
                << error_du_y[step] << std::endl;
  
  MPI_Finalize ();
  
  return 0;
}
