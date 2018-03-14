#include <tmesh_3d.h>
#include <simple_connectivity_3d_2trees.h>
#include <iostream>
int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;

  std::cout << "ciao" << std::endl << (int*)nullptr << std::endl
    << (int*)(tmsh.conn) << std::endl;
    
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  std::cout << "ciao" << std::endl << (int*)nullptr << std::endl
    << (int*)(tmsh.conn) << std::endl;
     
  return 0;
}
