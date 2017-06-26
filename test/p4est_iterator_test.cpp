#include <bim_timing.h>
#include <tmesh.h>

#include <mosfet_connectivity_2d.h>
#include <mosfet_doping_2d.h>

#include <vector>
#include <cassert>


static int
doping_driven_refinement (tmesh::quadrant_iterator quadrant)
{

  constexpr double L = 3.0e-6;
  constexpr double H = 1.0e-5;

  double maxy = 0, miny = 0, y = 0;
  double xcoord, ycoord;
  
  double top = std::numeric_limits<double>::lowest ();

  maxy = std::numeric_limits<double>::lowest ();
  miny = std::numeric_limits<double>::max ();

  for (int ii = 0; ii < 4; ++ii)
    {

      xcoord = (*quadrant).p(0, ii);
      ycoord = (*quadrant).p(1, ii);
      
      y = signedlog (doping (xcoord, ycoord, L, H));

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

  tmsh.set_refine_marker (doping_driven_refinement);

  recursive = 0;
  partforcoarsen = 0;

  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
    
  tmsh.vtk_export ("p4est_iterator_test");
  
  double xcoord = 0.0;
  double ycoord = 0.0;

  int ii = 0;
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      std::cout << ++ii;
      for (int ii = 0; ii < 4; ++ii)
        {
          xcoord = (*quadrant).p(0, ii);
          ycoord = (*quadrant).p(1, ii);
          std::cout << ", " << xcoord << ", " << ycoord;
        }
      std::cout << std::endl;
    }
  
  
  tmsh.vtk_export ("p4est_ref_test");

  MPI_Finalize ();
  return 0;

}

