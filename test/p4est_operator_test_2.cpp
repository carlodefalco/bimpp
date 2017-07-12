#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
#include <cassert>

static int
uniform_refinement (tmesh::quadrant_iterator quadrant)
{ return 1; }

int
main (int argc, char **argv)
{
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;
  
  MPI_Init (&argc, &argv);

  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_operator_test_2.octbin.gz");

  tmsh.read_connectivity ("p4est_operator_test_2.octbin.gz");
  
  // Uniform refinement.
  tmsh.set_refine_marker (uniform_refinement);
  
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_operator_test_2");
  
  // Assemble advection-diffusion matrix.
  sparse_matrix A;
  A.resize(tmsh.num_owned_nodes());
  
  double lambda = 25;
  std::vector<double> alpha(tmsh.num_local_elems (), 1);
  std::vector<double> psi(tmsh.num_local_nodes (), 0);
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          psi[quadrant->t(ii)] = lambda * (quadrant->p(0, ii) + quadrant->p(1, ii));
        }
    }
  
  bim2a_advection_diffusion (tmsh, alpha, psi, A);
  
  // Assemble right-hand side.
  std::vector<double> rhs(tmsh.num_local_nodes (), 0);
  
  std::vector<double> f(tmsh.num_local_elems (), 0);
  std::vector<double> g(tmsh.num_local_nodes (), 0);
  
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
  std::vector<double> vals;
  std::vector<int> irow, jcol;
  
  A.aij(vals, irow, jcol, 1);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  std::cout << "Solving linear system." << std::endl;
  
  mumps mumps_solver;
  
  if (rank == 0)
    mumps_solver.set_lhs_structure(A.rows(), irow, jcol);
  
  mumps_solver.analyze();
  
  if (rank == 0)
    mumps_solver.set_lhs_data(vals);
  
  mumps_solver.factorize();
  
  if (rank == 0)
    mumps_solver.set_rhs(rhs);
  
  mumps_solver.solve();
  
  // Compute error.
  double local_error = 0;
  double error = std::numeric_limits<double>::lowest ();
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          local_error = std::abs(u_ex(quadrant->p(0, ii), quadrant->p(1, ii)) - rhs[quadrant->t(ii)]);
          error = error < local_error ? local_error : error;
        }
    }
  
  std::cout << "L^inf-norm error: " << error << "." << std::endl;
  
  // Export solution.
  tmsh.octbin_export ("p4est_operator_test_2_output", rhs);
  
  MPI_Finalize ();
  
  return 0;
}
