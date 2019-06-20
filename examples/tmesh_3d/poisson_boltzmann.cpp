#include <bim_timing.h>
#include <tmesh_3d.h>
#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>

#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include <p8est.h>

constexpr double ll = -12;
constexpr double rr =  30;
constexpr double bb = -12;
constexpr double tt =  30;

constexpr p4est_topidx_t simple_conn_num_vertices = 8;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*3] = 
  {ll, ll, bb, rr, ll, bb, ll, rr, bb, rr, rr, bb,
   ll, ll, tt, rr, ll, tt, ll, rr, tt, rr, rr, tt};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*9] = 
  {1, 2, 3, 4, 5, 6, 7, 8, 1};


static constexpr double stern_size  = rr / 1.e-2 ;
static constexpr int maxlevel       =  10;
static constexpr int minlevel       =  2;
static constexpr double decay       = -3.0;
static char filename[255];

struct
ion
{
  std::array<double, 3> center;
  double radius;

  ion (const std::array<double, 3>& c,
       const double& r)
    : center (c), radius (r) {};
};

static std::vector<ion> ions;

double
levelsetfun (double x, double y, double z)
{
  double dist = 0.0; 
  for (const ion& i : ions)
    {
      
      dist += std::exp (decay * ((std::pow (x - i.center[0], 2) +
                                  std::pow (y - i.center[1], 2) +
                                  std::pow (z - i.center[2], 2)) /
                                 std::pow (i.radius, 2) - 1.0));            
      
    }
  return dist;
}


static int
uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
{ return 1; }

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

#include "crambina.h"
  
  TIC ();
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  TOC ("read connectivity");
  
  for (auto i = 0; i < minlevel; ++i)
    {

      TIC ();
      if (rank == 0) { tic (); }

      tmsh.set_refine_marker (uniform_refinement);
      recursive = 0; partforcoarsen = 1;
      tmsh.refine (recursive, partforcoarsen);

      TOC ("uniform refinement");
      
    }

  TIC ();
  sprintf (filename, "poisson_boltzmann_initial_mesh");
  tmsh.vtk_export (filename);
  TOC ("save initial mesh");


  for (int kk = 0; kk < (maxlevel - minlevel); ++kk)
    {

      {
        TIC();
        distributed_vector rcoeff (tmsh.num_owned_nodes ());
  
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant)
          {

            for (int ii = 0; ii < 8; ++ii)
              {
                if (! quadrant->is_hanging (ii))
                  rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                                           quadrant->p (1, ii),
                                                           quadrant->p (2, ii));
              }
          }
        bim3a_solution_with_ghosts (tmsh, rcoeff, replace_op);
        TOC ("compute coefficient");

  
        TIC(); 
        auto refinement = [&rcoeff]
          (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 0;
            double min = 100.0 * ions.size ();
            double max = 0.0;
            double tmp = 0.0;

            if (currentlevel >= maxlevel)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {

                    if (! q->is_hanging (ii))
                      tmp = rcoeff[q->gt (ii)];

                    if (tmp > max) max = tmp;
                    if (tmp < min) min = tmp;
     
                  }
                if (max >= 1.0 && min <= 1.0)
                  retval = maxlevel - currentlevel;
              }

            return (retval);
          };

        tmsh.set_refine_marker (refinement);
        recursive = 0; partforcoarsen = 1; balance = 0;
        tmsh.refine (recursive, partforcoarsen);      
        TOC ("refinement");
      }

      {
        TIC();
        distributed_vector rcoeff (tmsh.num_owned_nodes ());
  
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant)
          {

            for (int ii = 0; ii < 8; ++ii)
              {
                if (! quadrant->is_hanging (ii))
                  rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p(0, ii),
                                                           quadrant->p(1, ii),
                                                           quadrant->p(2, ii));
                else
                  for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                    rcoeff[quadrant->gparent (jj, ii)] += 0;
              }
          }
        bim3a_solution_with_ghosts (tmsh, rcoeff, replace_op);
        TOC ("compute coefficient");

        TIC ();
        sprintf (filename, "poisson_boltzmann_mesh_marked_%4.4d", kk);
        tmsh.octbin_export (filename, rcoeff);
        TOC ("save marker");

        TIC(); 
        auto coarsening = [&rcoeff]
          (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 0;
            double min = 100.0 * ions.size ();
            double max = 0.0;
            double tmp = 0.0;

            if (currentlevel <= minlevel)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {

                    if (! q->is_hanging (ii))
                      tmp = rcoeff[q->gt (ii)];

                    if (tmp > max) max = tmp;
                    if (tmp < min) min = tmp;
     
                  }
                if (max < 1.0 || min > 1.0)
                  retval = currentlevel - minlevel;
              }

            return (retval);
          };

        
        tmsh.set_coarsen_marker (coarsening);
        recursive = 0; partforcoarsen = 1;
        tmsh.coarsen (recursive, partforcoarsen);      
        TOC ("coarsening");
       
      }
    }

  {
    TIC();
    distributed_vector rcoeff (tmsh.num_owned_nodes ());
  
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
      {
            
        for (int ii = 0; ii < 8; ++ii)
          {
            if (! quadrant->is_hanging (ii))
              rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p(0, ii),
                                                       quadrant->p(1, ii),
                                                       quadrant->p(2, ii));
            else
              for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                rcoeff[quadrant->gparent (jj, ii)] += 0.;
          }
      }
    bim3a_solution_with_ghosts (tmsh, rcoeff, replace_op);
    TOC ("compute coefficient");

    TIC ();
    sprintf (filename, "poisson_boltzmann_mesh_marked_%4.4d", (maxlevel - minlevel));
    tmsh.octbin_export (filename, rcoeff);
    TOC ("save marker");
  }
  
  if (rank == 0) {print_timing_report();}
  
  MPI_Barrier (MPI_COMM_WORLD);
  
  MPI_Finalize ();
  return 0;

}

