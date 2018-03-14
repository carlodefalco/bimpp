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
uniform_refinement (tmesh::quadrant_iterator quadrant)
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
  
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 8; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  tmsh.set_refine_marker (top_refinement);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.set_refine_marker (top_refinement);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.set_refine_marker (right_refinement);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.set_refine_marker (right_refinement);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_operator_test_2");
  
  // Assemble advection-diffusion matrix.
  sparse_matrix A;
  A.resize(tmsh.num_global_nodes());
  
  double lambda = 25;
  std::vector<double> alpha(tmsh.num_local_quadrants (), 1);
  std::vector<double> psi(tmsh.num_global_nodes (), 0);
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
           if (! quadrant->is_hanging (ii))
            {
              psi[quadrant->gt(ii)] =
                lambda * (quadrant->p(0, ii) + quadrant->p(1, ii));
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
  func u_ex =
    [lambda] (double x, double y)
    { return (exp(lambda * x) - 1) / (exp(lambda) - 1) *
             (exp(lambda * y) - 1) / (exp(lambda) - 1); };
             
  dirichlet_bcs bcs;
  for (int i = 0; i < 4; ++i)
    bcs.push_back (std::make_tuple(0, i, u_ex));
  
  bim2a_dirichlet_bc (tmsh, bcs, A, rhs);
  
  // Solve problem.
  std::cout << "Solving linear system." << std::endl;
  
  mumps mumps_solver;
  
  std::vector<double> vals;
  std::vector<int> irow, jcol;
  
  A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
  
  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);
  
  // Reduce rhs (so that rank 0 has the actual rhs).
  std::vector<double> global_rhs(tmsh.num_global_nodes(), 0);
  MPI_Allreduce(rhs.data(), global_rhs.data(), rhs.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (rank == 0)
    mumps_solver.set_rhs (global_rhs);
  
  // Solve.
  mumps_solver.analyze ();
  mumps_solver.factorize ();
  mumps_solver.solve ();
  mumps_solver.cleanup ();
  
  // Export solution.
  MPI_Bcast(global_rhs.data(), global_rhs.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  tmsh.octbin_export ("p4est_operator_test_2_output", global_rhs);
  
  MPI_Finalize ();
  
  return 0;
}
