#include <bim_timing.h>
#include <mosfet_doping_3d.h>
#include <mosfet_connectivity_3d.h>
#include <tmesh_3d.h>

#include <vector>
#include <cassert>

static int
doping_driven_refinement (tmesh_3d::quadrant_iterator quadrant)
{

  constexpr double L = 3.0e-6;
  constexpr double H = 1.0e-5;
  constexpr double W = 1.0e-6;

  double maxy = 0, miny = 0, y = 0;
  double xcoord, ycoord, zcoord;
  
  double top = std::numeric_limits<double>::lowest ();

  maxy = std::numeric_limits<double>::lowest ();
  miny = std::numeric_limits<double>::max ();

  for (int ii = 0; ii < 8; ++ii)
    {

      xcoord = quadrant->p(0, ii);
      ycoord = quadrant->p(1, ii);
      zcoord = quadrant->p(2, ii);
      
      y = signedlog (doping (xcoord, ycoord, zcoord, L, H, W));

      maxy = maxy < y ? y : maxy;
      miny = miny > y ? y : miny;
      top  = top < ycoord ? ycoord : top;
    }

  double delta = maxy - miny;
  return ((top <= 0 && delta > .1) ? 1 : 0);
  
}


int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity3 ("p4est_ref_test_3d.octbin.gz");

  tmsh.read_connectivity ("p4est_ref_test_3d.octbin.gz");

  recursive = 0;
  partforcoarsen = 0;
 
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  for (int k = 0; k < 13 ; ++k)
    {
      tmsh.set_refine_marker (doping_driven_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("refinement and balancing"); }


  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  tmsh.vtk_export ("p4est_ref_test_3d");

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("IO"); }


  if (rank == 0)
    { print_timing_report (); }

  MPI_Finalize ();
  return 0;

}
