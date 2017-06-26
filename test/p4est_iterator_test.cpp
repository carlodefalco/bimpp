#include <bim_timing.h>
#include <mosfet_connectivity_2d.h>
#include <tmesh.h>

#include <vector>
#include <cassert>

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
    write_example_connectivity ("p4est_ref_test.octbin.gz");

  tmsh.read_connectivity ("p4est_ref_test.octbin.gz");

  double xcoord = 0.0;
  double ycoord = 0.0;

  for (auto quadrant = tmsh.begin_quadrant_sweep();
       quadrant != tmsh.end_quadrant_sweep();
       ++quadrant)
  {
    /*for (int ii = 0; ii < 4; ++ii)
    {
      xcoord = (*quadrant).p(0, ii);
      ycoord = (*quadrant).p(1, ii);
    }*/
  }
  
  tmsh.vtk_export ("p4est_ref_test");

  MPI_Finalize ();
  return 0;

}

