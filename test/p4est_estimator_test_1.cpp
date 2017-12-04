#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
#include <cassert>

static int
top_refinement (tmesh::quadrant_iterator quadrant)
{
  double ycoord;
  double bottom = std::numeric_limits<double>::max ();
  for (int ii = 0; ii < 4; ++ii)
  {
    ycoord = quadrant->p(1, ii);
    bottom = bottom > ycoord ? ycoord : bottom;
  }
  return ((bottom >= 0.9) ? 1 : 0);
}

static int
right_refinement (tmesh::quadrant_iterator quadrant)
{
  double xcoord;
  double left = std::numeric_limits<double>::max ();
  for (int ii = 0; ii < 4; ++ii)
  {
    xcoord = quadrant->p(0, ii);
    left = left > xcoord ? xcoord : left;
  }
  return ((left >= 0.9) ? 1 : 0);
}

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
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  tmsh.set_refine_marker (uniform_refinement);
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 2; ++cycle)
    tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_estimator_test_1");
  
  for (int adapt = 0; adapt < 7; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
      // Assemble matrix.
      sparse_matrix A, M;
      A.resize(tmsh.num_global_nodes());
      M.resize(tmsh.num_global_nodes());
      
      double epsilon = 1e-3;
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
      tmsh.octbin_export ((std::string("p4est_estimator_test_1_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::cout << " Done." << std::endl;
      
      // Compute reconstructed gradient.
      std::cout << "Computing reconstructed gradient and estimator.";
      
      gradient du = bim2c_quadtree_pde_recovered_gradient(tmsh, global_rhs);
      
      auto refine_fun = [& du, & global_rhs, &tmsh] (tmesh::quadrant_iterator q)
        { return zz_marker_grad (q, du, global_rhs,
                                 1e-2 / std::sqrt(tmsh.num_global_nodes())); };
      
      tmsh.octbin_export ((std::string("p4est_estimator_test_1_du_x_")
                           + std::to_string(adapt)).c_str(), du.first);
      tmsh.octbin_export ((std::string("p4est_estimator_test_1_du_y_")
                           + std::to_string(adapt)).c_str(), du.first);
      
      // Refine according to refine_fun.
      tmsh.set_refine_marker (refine_fun);
      tmsh.refine (recursive, partforcoarsen);
      
      tmsh.vtk_export ((std::string("p4est_estimator_test_1_refined_")
                        + std::to_string(adapt)).c_str());
      std::cout << " Done." << std::endl;
    }
  
  MPI_Finalize ();
  
  return 0;
}
