#include <mosfet_connectivity_2d.h>
#include <tmesh.h>

#include <vector>
#include <cassert>


int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh, tmsh2;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_load_save_test2.octbin.gz");

  tmsh.read_connectivity ("p4est_load_save_test2.octbin.gz");  
  
  std::vector<double> prova(tmsh.num_global_nodes(), 0);
  tmsh.octbin_export ("p4est_load_save_test2out.octbin.gz",prova);

  MPI_Finalize ();
  return 0;

}
