#include <bim_sparse.h>
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
    write_example_connectivity ("p4est_operator_test.octbin.gz");

  tmsh.read_connectivity ("p4est_operator_test.octbin.gz");
  
  // Uniform refinement.
  tmsh.set_refine_marker (uniform_refinement);
  
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_operator_test");
  
  // Assemble advection-diffusion matrix.
  sparse_matrix A;
  A.resize(tmsh.num_owned_nodes());
  
  std::vector<double> alpha(tmsh.num_local_elems (), 1);
  std::vector<double> psi(tmsh.num_local_nodes (), 0);
  
  // alpha = x; psi(x, y) = x * y;
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      alpha.push_back(1);// 0.5 * (quadrant->p(0, 0) + quadrant->p(0, 1)) );
      
      for (int ii = 0; ii < 4; ++ii)
        {
          psi[quadrant->t(ii)] = 0;
        }
    };
  
  bim2a_advection_diffusion (tmsh, alpha, psi, A);
  
  std::cout << A << std::endl;
  
  MPI_Finalize ();
  
  return 0;
}
