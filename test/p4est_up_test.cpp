#include <bim_timing.h>
#include <simple_connectivity_2d_2trees.h>
#include <tmesh.h>

#include <vector>
#include <cassert>

static int
fake_refinement (tmesh::quadrant_iterator quadrant)
{
  static int pippo = 0;
  std::cout << "call #" << pippo++ << " level #"
            << (int) quadrant->the_quadrant->level << " tree #"
            << (int) quadrant->get_tree_idx () << " forest_quad_idx #"
            << (int) quadrant->get_forest_quad_idx () << " tree_quad_idx #"
            << (int) quadrant->get_tree_quad_idx () 
            << std::endl;
  return 1;
}

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int         recursive, partforcoarsen, balance;
  MPI_Comm    mpicomm = MPI_COMM_WORLD;  
  int         rank, size;
  tmesh       tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  tmsh.set_refine_marker (fake_refinement);

  recursive = 0;
  partforcoarsen = 0;
 
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  for (int k = 0; k < 5 ; ++k)
    {
      std::cout << "refinement step #" << k << std::endl;
      tmsh.refine (recursive, partforcoarsen);
    }

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("refinement and balancing"); }


  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  tmsh.vtk_export ("p4est_ref_test");

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("IO"); }


  if (rank == 0)
    { print_timing_report (); }

  MPI_Finalize ();
  return 0;

}

