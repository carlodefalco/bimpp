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
    write_example_connectivity ("p4est_iterator_test.octbin.gz");

  tmsh.read_connectivity ("p4est_iterator_test.octbin.gz");

  tmsh.set_refine_marker (uniform_refinement);
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);

  std::cout << "first refinement step"
            << std::endl;
  std::cout << "num_local_nodes = "
            << tmsh.num_local_nodes ()
            << std::endl;

  std::cout << "num_owned_nodes = "
            << tmsh.num_owned_nodes ()
            << std::endl;

  tmsh.set_refine_marker (bottom_refinement);
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);

  std::cout << "second refinement step"
            << std::endl;
  std::cout << "num_local_nodes = "
            << tmsh.num_local_nodes ()
            << std::endl;

  std::cout << "num_owned_nodes = "
            << tmsh.num_owned_nodes ()
            << std::endl;

  tmsh.set_refine_marker (left_refinement);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  
  std::cout << "third refinement step"
            << std::endl;
  std::cout << "num_local_nodes = "
            << tmsh.num_local_nodes ()
            << std::endl;

  std::cout << "num_owned_nodes = "
            << tmsh.num_owned_nodes ()
            << std::endl;
    
  int ii = 0;
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      std::cout << ++ii;
      for (int ii = 0; ii < 4; ++ii)
        {
          std::cout << ", ";
          if (quadrant->is_hanging (ii))
            std::cout << "*";
          else
            std::cout << quadrant->t(ii);
        }
      std::cout << std::endl;
    }
  
  
  tmsh.vtk_export ("p4est_iterator_test");

  MPI_Finalize ();
  return 0;

}

