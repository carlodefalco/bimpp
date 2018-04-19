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
static auto zero = [] (double, double) {return 0;};

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
  
  tmsh.vtk_export ("p4est_adr_test_2_isogeo");
  
  std::vector<tmesh::idx_t> nnodes;
  std::vector<double> h_step;
  
  double delta1 = 0.01;
  double delta2 = 0.01;
  
  for (int adapt = 0; adapt < refine_steps; ++adapt)
    {
      std::cout << "*** Step " << adapt << " ***" << std::endl;
      
      // Assemble matrix.
      sparse_matrix A;
      A.resize(tmsh.num_global_nodes());
      
      double epsilon = 1e-6;
      double theta = M_PI / 4;
      std::vector<double> alpha(tmsh.num_local_quadrants (), epsilon);
      std::vector<double> psi(tmsh.num_global_nodes (), 0);
      
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
                }
            }
        }
      
      // Reduce coefficients.
      std::vector<double> global_psi(tmsh.num_global_nodes(), 0);
      MPI_Allreduce(psi.data(), global_psi.data(), psi.size(),
                    MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
      bim2a_advection_diffusion (tmsh, alpha, global_psi, A);
      
      // Assemble right-hand side.
      std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      
      std::vector<double> f(tmsh.num_local_quadrants (), 0);
      std::vector<double> g(tmsh.num_global_nodes (), 0);
      
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
      tmsh.octbin_export ((std::string("p4est_adr_test_2_isogeo_u_")
                           + std::to_string(adapt)).c_str(), global_rhs);
      
      std::cout << " Done." << std::endl;
      
      // Compute gradient norm and h.
      std::cout << "Computing gradient norm.";
      
      // Compute h.
      double hx = 0, hy = 0,
             h = std::numeric_limits<double>::max (),
             global_h = 0;
      
      double grad_sum = 0, global_grad = 0;
      
      std::vector<double> grad(tmsh.num_local_quadrants ());
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          hx = quadrant->p(0, 1) - quadrant->p(0, 0);
          hy = quadrant->p(1, 2) - quadrant->p(1, 0);
          
          h = std::min(h, std::sqrt(hx*hx + hy*hy));
          
          grad[quadrant->get_forest_quad_idx ()]
            = semih1_error (quadrant, zero, zero, global_rhs) / (hx * hy);
          
          grad_sum += grad[quadrant->get_forest_quad_idx ()];
        }
      
      MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);
      MPI_Allreduce(&grad_sum, &global_grad, 1, MPI_DOUBLE, MPI_SUM, mpicomm);
      
      nnodes.push_back (tmsh.num_global_nodes ());
      h_step.push_back (global_h);
      
      std::cout << " Done." << std::endl;
      
      auto refine_fun = [& delta1, & grad, & global_grad, & tmsh] (tmesh::quadrant_iterator q)
        {
          return grad[q->get_forest_quad_idx ()] >=
                   delta1 * global_grad / tmsh.num_global_quadrants ();
        };
      
      auto coarsen_fun = [& delta2, & grad, & global_grad, & tmsh] (tmesh::quadrant_iterator q)
        {
          return grad[q->get_forest_quad_idx ()] <=
                   delta2 * global_grad / tmsh.num_global_quadrants ();
        };
      
      if (tmsh.num_global_nodes () >= 1e6)
        break;
      
      // Coarsen and refine.
      tmsh.set_coarsen_marker (coarsen_fun);
      tmsh.set_refine_marker (refine_fun);
      
      tmsh.coarsen (recursive, partforcoarsen, 0);
      tmsh.refine (recursive, partforcoarsen);
      
      tmsh.vtk_export ((std::string("p4est_adr_test_2_isogeo_newmesh_")
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
