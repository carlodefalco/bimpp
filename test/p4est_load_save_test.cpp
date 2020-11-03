#include <bim_timing.h>
#include <mosfet_doping_2d.h>
#include <mosfet_connectivity_2d.h>
#include <tmesh.h>

#include <vector>

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh, tmsh2;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_load_save_test.octbin");

  tmsh.read_connectivity ("p4est_load_save_test.octbin.gz");  
  tmsh.save ("p4est_load_save_test.p4est");  
  tmsh2.load ("p4est_load_save_test.p4est");
  tmsh.vtk_export ("p4est_load_save_test");
  tmsh2.vtk_export ("p4est_load_save_test_2");

  MPI_Finalize ();
  return 0;

}

