#include <bim_timing.h>
#include <mosfet_doping_2d.h>
#include <mosfet_connectivity_2d.h>
#include <tmesh.h>

#include <vector>
#include <cassert>

/*
static int
refine_fn (p4est_t * p4est, p4est_topidx_t tt,
           p4est_quadrant_t * quadrant)
{
  p4est_quadrant_t node;
  int i;
  double vxyz[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (i = 0; i < P4EST_CHILDREN; ++i)
    {
      p4est_quadrant_corner_node (quadrant, i, &node);
      p4est_qcoord_to_vertex (p4est->connectivity, tt, node.x, node.y,
                              &(vxyz[3 * i]));
    }

  return (doping_driven_refinement (vxyz));

}
*/

static int
refine_fn (tmesh::quadrant_t &quadrant)
{
  return (doping_driven_refinement ([&quadrant] (tmesh::idx_t i, tmesh::idx_t j)
                                    {return quadrant.p (i, j); }));
}

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

  recursive = 0;
  partforcoarsen = 0;

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  for (int k = 0; k < 13 ; ++k)
    {
      p4est_refine (tmsh.p4est, recursive, refine_fn, NULL);
      p4est_balance (tmsh.p4est, P4EST_CONNECT_FACE, NULL);
      p4est_partition (tmsh.p4est, partforcoarsen, NULL);
    }

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("refinement and balancing"); }


  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  p4est_vtk_write_file (tmsh.p4est, NULL, "p4est_ref_test");

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("IO"); }


  if (rank == 0)
    { print_timing_report (); }

  MPI_Finalize ();
  return 0;

}

