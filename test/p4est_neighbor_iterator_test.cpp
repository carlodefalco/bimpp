#include <bim_timing.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
#include <cassert>


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

  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;
  double                xcoord = 0.0;
  double                ycoord = 0.0;

  MPI_Init (&argc, &argv);

  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_neighbor_iterator_test.octbin.gz");

  tmsh.read_connectivity ("p4est_neighbor_iterator_test.octbin.gz");

  tmsh.set_refine_marker (uniform_refinement);
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);

  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (auto neighbor = quadrant->begin_neighbor_sweep ();
           neighbor != quadrant->end_neighbor_sweep ();
           ++neighbor)
        {
          std::cout << "Element " << quadrant->get_global_quad_idx()
                    << ", neighbor " << neighbor->get_global_quad_idx()
                    << ", vertices: " << neighbor->t(0) << ", " << neighbor->t(1)
                    << ", " << neighbor->t(2) << ", " << neighbor->t(3)
                    << std::endl;

        }
      
      std::cout << std::endl;
    }
  
  tmsh.vtk_export ("p4est_neighbor_iterator_test");

  MPI_Finalize ();
  return 0;

}

