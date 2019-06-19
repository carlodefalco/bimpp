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

constexpr double ll = -.3e-9;
constexpr double rr =  .3e-9;
constexpr double bb = -.3e-9;
constexpr double tt =  .3e-9;

constexpr p4est_topidx_t simple_conn_num_vertices = 8;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*3] = 
  {ll, ll, bb, rr, ll, bb, ll, rr, bb, rr, rr, bb,
   ll, ll, tt, rr, ll, tt, ll, rr, tt, rr, rr, tt};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*9] = 
  {1, 2, 3, 4, 5, 6, 7, 8, 1};


static constexpr double stern_size  = rr / 1.e-2 ;
static constexpr int maxlevel       = 10;
static constexpr int minlevel       = 7;
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



static double
refinement (tmesh_3d::quadrant_iterator quadrant)
{

  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double xcoord, ycoord, zcoord;
  int retval = 0;
  marker m = exterior, m0 = exterior;
  
  for (int ii = 0; ii < 8; ++ii)
    {

      m = mark_region (quadrant->p(0, ii),
                       quadrant->p(1, ii),
                       quadrant->p(2, ii));
      if (ii == 0)
        m0 = m;
      
      if ((m == stern)
          || (m != m0))
        {
          retval = maxlevel - currentlevel;
          break;
        }
    }

  if (currentlevel >= maxlevel)
    retval = 0;
      
  return (retval);
}

static int
coarsening (tmesh_3d::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double xcoord, ycoord, zcoord;
  int retval = currentlevel - minlevel;
  marker m = exterior;
  
  for (int ii = 0; ii < 8; ++ii)
    {

      m = mark_region (quadrant->p(0, ii),
                       quadrant->p(1, ii),
                       quadrant->p(2, ii));


      if ((m != exterior)
          && (m != interior))
        {
          retval = 0;
          break;
        }
    }

  if (currentlevel <= minlevel)
    retval = 0;
  
  return (retval);
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


  ions.push_back (ion ({ 1.390000000000000e-10,  0.000000000000000e+00, 0.000000000000000e+00}, 7.000000000000001e-11));
  ions.push_back (ion ({ 6.950000000000002e-11,  1.203775311260370e-10, 0.000000000000000e+00}, 7.000000000000001e-11));
  ions.push_back (ion ({-6.949999999999998e-11,  1.203775311260370e-10, 0.000000000000000e+00}, 7.000000000000001e-11));
  ions.push_back (ion ({-1.390000000000000e-10,  0.000000000000000e+00, 0.000000000000000e+00}, 7.000000000000001e-11));
  ions.push_back (ion ({-6.950000000000007e-11, -1.203775311260369e-10, 0.000000000000000e+00}, 7.000000000000001e-11));
  ions.push_back (ion ({ 6.949999999999990e-11, -1.203775311260370e-10, 0.000000000000000e+00}, 7.000000000000001e-11));
  ions.push_back (ion ({ 2.480000000000000e-10,  0.000000000000000e+00, 0.000000000000000e+00}, 2.500000000000001e-11));
  ions.push_back (ion ({ 1.240000000000000e-10,  2.147743001385408e-10, 0.000000000000000e+00}, 2.500000000000001e-11));
  ions.push_back (ion ({-1.240000000000000e-10,  2.147743001385408e-10, 0.000000000000000e+00}, 2.500000000000001e-11));
  ions.push_back (ion ({-2.480000000000000e-10,  0.000000000000000e+00, 0.000000000000000e+00}, 2.500000000000001e-11));
  ions.push_back (ion ({-1.240000000000001e-10, -2.147743001385407e-10, 0.000000000000000e+00}, 2.500000000000001e-11));
  ions.push_back (ion ({ 1.239999999999998e-10, -2.147743001385409e-10, 0.000000000000000e+00}, 2.500000000000001e-11));
  
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
  TOC ("i/o");
  
  TIC ();
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

  
  TIC (); 
  auto estimator = [&rcoeff]
        (tmesh_3d::quadrant_iterator q)
        {
          if (rho2(q->centroid(0), q->centroid(1), q->centroid(2)) > (R*R))
            return estimator_sol (q, u_star0, result);
          else
            return estimator_sol (q, u_star1, result);
        };

  TOC ("refinement"); 

  /*
  TIC (); 
  recursive = 1;  partforcoarsen = 1;  
  tmsh.set_coarsen_marker (coarsening);
  tmsh.coarsen (recursive, partforcoarsen);
  TOC ("coarsening"); 
  
  TIC (); 
  sprintf (filename, "poisson_boltzmann_adapted_mesh");
  tmsh.vtk_export (filename);
  TOC ("i/o");
  */

  TIC (); 
  tmsh.octbin_export ("poisson_boltzmann_mesh_marked_0000", rcoeff);
  TOC ("i/o");
  
  if (rank == 0) {print_timing_report();}
  
  MPI_Barrier (MPI_COMM_WORLD);
  
  MPI_Finalize ();
  return 0;

}

