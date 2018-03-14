#include <mosfet_connectivity_3d.h>
#include <tmesh_3d.h>

#include <vector>
#include <cassert>


int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh, tmsh2;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity3 ("p8est_load_save_test2.octbin.gz");

  tmsh.read_connectivity ("p8est_load_save_test2.octbin.gz");  
  
  std::vector<double> prova(tmsh.num_global_nodes(), 0);
  tmsh.save ("p8est_load_save_test2.p8est");  
  tmsh2.load ("p8est_load_save_test2.p8est");
  tmsh.vtk_export ("p8est_load_save_test2");
  tmsh2.vtk_export ("p8est_load_save_test2_2");
  //tmsh.octbin_export ("p8est_load_save_test2",prova);

  MPI_Finalize ();
  return 0;

}
