#include <bim_timing.h>
#include <mosfet_doping_2d.h>
#include <mosfet_connectivity_2d.h>
#include <tmesh.h>

#include <vector>
#include <cassert>
#include <cstdio>

constexpr double L = 3.0e-6;
constexpr double H = 1.0e-5;

static int
doping_driven_refinement (tmesh::quadrant_iterator quadrant,
                          const std::vector<double> &y,
                          const std::vector<double> &ycoord)
{
  double maxy = std::numeric_limits<double>::lowest ();
  double miny = std::numeric_limits<double>::max ();
  double top  = std::numeric_limits<double>::lowest ();

  for (int ii = 0; ii < 4; ++ii)
    {
      double tmp = 0;
      
      if (! quadrant->is_hanging(ii))
        tmp = y[quadrant->t(ii)];
      else
        tmp = (y[quadrant->parent(0, ii)] +
               y[quadrant->parent(1, ii)]) / 2.;
          
      maxy = maxy < tmp ? tmp : maxy;
      miny = miny > tmp ? tmp : miny;

      if (! quadrant->is_hanging(ii))
        tmp = ycoord[quadrant->t(ii)];
      else
        tmp = (ycoord[quadrant->parent(0, ii)] +
               ycoord[quadrant->parent(1, ii)]) / 2.;

      top  = top < tmp ? tmp : top;   
    }

  double delta = maxy - miny;
  return ((top <= 0 && delta > .1) ? 1 : 0);
  
}

static int
coarsen_right_half (std::function<tmesh::quadrant_iterator (tmesh::idx_t)> next,
                    const std::vector<double> &xcoord)
{
  double left = L;
  tmesh::quadrant_iterator qi;
  for (tmesh::idx_t ii = 0; ii < 4; ++ii)
    {
      qi = next (ii);
      auto ud = (tmesh_qdata_t *) (qi->the_quadrant->p.user_data);
      for (int jj = 0; jj < 4; ++jj)
        {
          double tmp = 0;
          
          if (! ud->is_hanging[ii])
            tmp = xcoord[ud->t[ii]];
          else
            tmp = (xcoord[ud->t[ii]] +
                   xcoord[ud->t[ii+4]]) / 2.;

          left  = left > tmp ? tmp : left;   
        }     
    }
  return (left > L/2.0 ? 1 : 0);
}

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int       recursive, partforcoarsen, balance;
  MPI_Comm  mpicomm = MPI_COMM_WORLD;  
  int       rank, size;
  tmesh     tmsh;
  char      filename[128];
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_ref_lnodes_test.octbin.gz");

  tmsh.read_connectivity ("p4est_ref_lnodes_test.octbin.gz");

  recursive = 0;
  partforcoarsen = 1;
 

  for (int k = 0; k < 22 ; ++k)
    {
      std::vector<double>
        y(tmsh.num_local_nodes ()),
        ycoord(tmsh.num_local_nodes ()),
        xcoord(tmsh.num_local_nodes ());


      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          for (int jj = 0; jj < 4; ++jj)
            if (! quadrant->is_hanging (jj))
              {
                ycoord[quadrant->t(jj)] = quadrant->p(1, jj);
                xcoord[quadrant->t(jj)] = quadrant->p(0, jj);
                y[quadrant->t(jj)] =
                  signedlog (doping (xcoord[quadrant->t(jj)],
                                     ycoord[quadrant->t(jj)],
                                     L, H));                
              }
        }


      if (k < 11)
        {

          MPI_Barrier (MPI_COMM_WORLD);
          if (rank == 0)
            { tic (); }
          
          auto
            refmark = [&y, &ycoord]
            (tmesh::quadrant_iterator quadrant) -> int
            {
              return doping_driven_refinement (quadrant, y, ycoord);
            };
      
          tmsh.set_refine_marker (refmark);
          tmsh.update ();
          tmsh.refine (recursive, partforcoarsen);

            MPI_Barrier (MPI_COMM_WORLD);
            if (rank == 0)
              { toc ("refinement and balancing"); }
        }
      else
        {

          MPI_Barrier (MPI_COMM_WORLD);
          if (rank == 0)
            { tic (); }
          
          auto
            coamark = [&xcoord]
            (std::function<tmesh::quadrant_iterator (tmesh::idx_t)> next) -> int
            {
              return coarsen_right_half (next, xcoord);
            };
      
          tmsh.set_coarsen_marker (coamark);
          tmsh.update ();
          tmsh.coarsen (recursive, partforcoarsen);

           MPI_Barrier (MPI_COMM_WORLD);
            if (rank == 0)
              { toc ("coarsening"); }
        }
      
    

      MPI_Barrier (MPI_COMM_WORLD);
      if (rank == 0)
        { tic (); }

      sprintf (filename, "p4est_ref_lnodes_test_%3.3d", k);
      tmsh.vtk_export (filename);

      MPI_Barrier (MPI_COMM_WORLD);
      if (rank == 0)
        { toc ("IO"); }

    }

  if (rank == 0)
    { print_timing_report (); }

  MPI_Finalize ();
  return 0;

}
