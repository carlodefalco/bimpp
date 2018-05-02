#include <iostream>
#include <cassert>
#include <vector>
#include <limits>

#include <bim_timing.h>
#include <tmesh.h>
#include <simple_connectivity_2d.h>

char filename[255];

static int
bottom_refinement (tmesh::quadrant_iterator quadrant)
{
  double ycoord;
  double top = std::numeric_limits<double>::lowest ();
  for (int ii = 0; ii < 4; ++ii)
    {
      ycoord = quadrant->p(1, ii);
      top  = top < ycoord ? ycoord : top;
    }
  return ((top <= 0.5) ? 1 : 0);
}

static int
left_refinement (tmesh::quadrant_iterator quadrant)
{
  double xcoord;
  double right = std::numeric_limits<double>::lowest ();
  for (int ii = 0; ii < 4; ++ii)
    {
      xcoord = quadrant->p(0, ii);
      right  = right < xcoord ? xcoord : right;
    }
  return ((right <= 0.5) ? 1 : 0);
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
  double                xcoord = 0.0;
  double                ycoord = 0.0;

  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  recursive = 0; partforcoarsen = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.set_refine_marker (bottom_refinement);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.set_refine_marker (bottom_refinement);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.set_refine_marker (left_refinement);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.set_refine_marker (left_refinement);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_quad_index_test");
  
  int iel = 0;

  for (int irank = 0; irank < size; ++irank)
    {
      if (irank == rank)
        {
          std::cout << "rank " << rank << ": " << std::endl;
          for (auto quadrant = tmsh.begin_quadrant_sweep ();
               quadrant != tmsh.end_quadrant_sweep ();
               ++quadrant)
            {
              std::cout << "\tquad "
                        << quadrant->get_global_quad_idx ()
                        << std::endl;
            }
        }
      MPI_Barrier (mpicomm);
    }
        
  MPI_Finalize ();
  return 0;
}
