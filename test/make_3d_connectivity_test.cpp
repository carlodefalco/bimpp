#include <bim_timing.h>
#include <tmesh_3d.h>

#include <vector>

static int
fake_refinement (tmesh_3d::quadrant_iterator quadrant)
{
  int         rank, size;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  static int pippo = 0;
  std::cout << "rank " << rank << std::endl;
  std::cout << "call #" << pippo++ << " level #"
            << (int) quadrant->the_quadrant->level << " tree #"
            << (int) quadrant->get_tree_idx () << " forest_quad_idx #"
            << (int) quadrant->get_forest_quad_idx () << " tree_quad_idx #"
            << (int) quadrant->get_tree_quad_idx () << " vertices (local numbering) "
            << (int) quadrant->t (0) << " " << (int) quadrant->t (1) << " "
            << (int) quadrant->t (2) << " " << (int) quadrant->t (3) << " "
            << (int) quadrant->t (4) << " " << (int) quadrant->t (5) << " "
            << (int) quadrant->t (6) << " " << (int) quadrant->t (7) << " "
            << " vertices (global numbering) "
            << (int) quadrant->gt (0) << " " << (int) quadrant->gt (1) << " "
            << (int) quadrant->gt (2) << " " << (int) quadrant->gt (3) << " "
            << (int) quadrant->gt (4) << " " << (int) quadrant->gt (5) << " "
            << (int) quadrant->gt (6) << " " << (int) quadrant->gt (7) << " "
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
  tmesh_3d    tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  const p4est_topidx_t num_trees[3] = {2, 3, 4};
  const double step[3] = {1., 1., 1.};
  double *p{nullptr};
  p4est_topidx_t num_vertices{0}, total_num_trees{0};
  p4est_topidx_t *t{nullptr};
  std::vector<std::pair<p4est_topidx_t, p4est_topidx_t>> bcells;
  make_connectivity_3d (num_trees, step, p, num_vertices, t, total_num_trees, bcells);
  tmsh.read_connectivity (p, num_vertices, t, total_num_trees);
  for (auto const & pippo : bcells) {
    std::cout << "bcells : " << pippo.first << " " << pippo.second << std::endl;
  }
  
  recursive = 0;
  partforcoarsen = 0;
 
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  for (int k = 0; k < 3 ; ++k)
    {
      std::cout << "refinement step #" << k << std::endl;
      tmsh.set_refine_marker (fake_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }


  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("refinement and balancing"); }


  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  tmsh.vtk_export ("make_3d_connectivity_test");

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("IO"); }


  if (rank == 0)
    { print_timing_report (); }

  MPI_Finalize ();
  return 0;

}
