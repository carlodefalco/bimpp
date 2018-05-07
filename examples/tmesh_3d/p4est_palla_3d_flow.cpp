#include <bim_timing.h>
#include <tmesh_3d.h>
#include <quad_operators_3d.h>
#include <mumps_class.h>

#include <array>
#include <cassert>
#include <cstdio>

#include "simple_connectivity_3d.h"

static std::array<double, 3> x0 = {.5, .5, .5};
static std::array<double, 3> x = x0;
static std::array<double, 3> v0 = {0.2, 0.0, 0.8};
static std::array<double, 3> v = v0;
static std::array<double, 3> L = {1., 1., 1.};
static std::array<double, 3> g = {0., -9.81, 0.};
static constexpr double r = .0625;
static constexpr double dt = .01;
static const int maxlevel =  8;
static const int minlevel =  3;
static char filename[255] = "\0";

static int
uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
{ return 1; }

static int
refinement (tmesh_3d::quadrant_iterator quadrant)
{
  int currentlevel = static_cast<int> (quadrant->the_quadrant->level);
  double xcoord, ycoord, zcoord, dist = .0;
  int retval = 0;
  for (int ii = 0; ii < 8; ++ii)
    {
      xcoord = quadrant->p(0, ii);
      ycoord = quadrant->p(1, ii);
      zcoord = quadrant->p(2, ii);
      
      dist = std::sqrt (std::pow (xcoord - x[0], 2) +
                        std::pow (ycoord - x[1], 2) +
                        std::pow (zcoord - x[2], 2));


      if (dist > .9 * r && dist < 1.1 * r)
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
  double xcoord, ycoord, zcoord, dist = .0;
  int retval = currentlevel - minlevel;
  for (int ii = 0; ii < 8; ++ii)
    {
      xcoord = quadrant->p(0, ii);
      ycoord = quadrant->p(1, ii);
      zcoord = quadrant->p(2, ii);
      
      dist = std::sqrt (std::pow (xcoord - x[0], 2) +
                        std::pow (ycoord - x[1], 2) +
                        std::pow (zcoord - x[2], 2));


      if (dist > .9 * r && dist < 1.1 * r)
        {
          retval = 0;
          break;
        }
    }

  if (currentlevel <= minlevel)
    retval = 0;
  
  return (retval);
}

static void
palla_helper (double _x0, double _v0, double _L,
              double _g, double _dt, double _r,
              double &_x, double &_v)
{
  _v = _v0 + _dt * _g;
  _x = _x0 + _dt * _v;
  bool collision = ((_x + _r) > _L) || ((_x - _r) < 0);
  if (collision)
    {
      _v *= -1.;
      double dxupright = ((_L - (_x + _r)) -
                          std::abs (_L - (_x + _r)));
      _x += dxupright;
      double dxdownleft = ((0. - (_x - _r)) +
                           std::abs (0 - (_x - _r)));
      _x += dxdownleft;
    }
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

  
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { tic (); }
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  if (rank == 0) { toc ("read connectivity"); }
  MPI_Barrier (MPI_COMM_WORLD);

  for (auto i = 0; i < minlevel; ++i)
    {

      MPI_Barrier (MPI_COMM_WORLD);
      if (rank == 0) { tic (); }

      tmsh.set_refine_marker (uniform_refinement);
      recursive = 0; partforcoarsen = 1;
      tmsh.refine (recursive, partforcoarsen);

      if (rank == 0) { toc ("uniform refinement"); }
      MPI_Barrier (MPI_COMM_WORLD);

      sprintf (filename, "palla_initial_mesh_%5.5d", i);
      tmsh.vtk_export (filename);
    }

  int iframe = 0;
  for (double t = .0; t < 10; t += dt)
    {

      palla_helper (x0[0], v0[0], L[0], g[0], dt, r, x[0], v[0]);
      palla_helper (x0[1], v0[1], L[1], g[1], dt, r, x[1], v[1]);
      palla_helper (x0[2], v0[2], L[2], g[2], dt, r, x[2], v[2]);
      
      x0 = x; v0 = v;

      // Assemble differential problem.
      double kG = 1e6;
      double kS = 1;
      
      std::vector<double> alpha(tmsh.num_local_quadrants (), kG);
      std::vector<double> psi(tmsh.num_global_nodes (), 0);
      
      std::vector<double> f(tmsh.num_local_quadrants (), 0);
      std::vector<double> g(tmsh.num_global_nodes (), 0);
      
      double xx = 0, yy = 0, zz = 0;
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
	  if (std::pow(quadrant->centroid(0) - x[0], 2) +
	      std::pow(quadrant->centroid(1) - x[1], 2) +
	      std::pow(quadrant->centroid(2) - x[2], 2) <=
	      std::pow(r, 2))
	    alpha[quadrant->get_forest_quad_idx ()] = kS;
	}

      sparse_matrix A;
      A.resize(tmsh.num_global_nodes());
      
      bim3a_advection_diffusion (tmsh, alpha, psi, A);
      
      std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      bim3a_rhs (tmsh, f, g, rhs);

      // Set boundary conditions.
      dirichlet_bcs3 bcs;
      bcs.push_back (std::make_tuple(0, 0, [] (double x, double y, double z) {return 1;}));
      bcs.push_back (std::make_tuple(0, 1, [] (double x, double y, double z) {return 0;}));
      
      bim3a_dirichlet_bc (tmsh, bcs, A, rhs);
      
      // Solve problem.
      std::cout << "Solving linear system.";
      
      mumps mumps_solver;
      
      std::vector<double> vals;
      std::vector<int> irow, jcol;
      
      A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
      mumps_solver.set_distributed_lhs_data (vals);
      
      // Reduce rhs (so that rank 0 has the actual rhs).
      std::vector<double> global_rhs(tmsh.num_global_nodes(), 0);
      MPI_Reduce(rhs.data(), global_rhs.data(), rhs.size(),
                 MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      
      if (rank == 0)
        mumps_solver.set_rhs (global_rhs);
      
      // Solve.
      mumps_solver.analyze ();
      mumps_solver.factorize ();
      mumps_solver.solve ();
      mumps_solver.cleanup ();
      
      // Export solution.
      MPI_Bcast(global_rhs.data(), global_rhs.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      tmsh.octbin_export ((std::string("p4est_palla_3d_flow_u_")
                           + std::to_string(t)).c_str(), global_rhs);
      
      std::cout << " Done." << std::endl;

      // Export level.
      std::vector<double> level (tmsh.num_local_quadrants ());
      
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
	{
	  level[quadrant->get_forest_quad_idx ()] = quadrant->the_quadrant->level;
	}
      
      tmsh.octbin_export_quadrant ((std::string("p4est_palla_3d_flow_level_")
				    + std::to_string(t)).c_str(), level);
      
      // Refine + coarsen.
      MPI_Barrier (MPI_COMM_WORLD); 
      if (rank == 0) { tic (); }
      recursive = 1;  partforcoarsen = 1;  
      tmsh.set_refine_marker (refinement);
      tmsh.refine (recursive, partforcoarsen);
      if (rank == 0) { toc ("refinement"); }

            
      MPI_Barrier (MPI_COMM_WORLD); 
      if (rank == 0) { tic (); }
      recursive = 1;  partforcoarsen = 1;  
      tmsh.set_coarsen_marker (coarsening);
      tmsh.coarsen (recursive, partforcoarsen);
      if (rank == 0) { toc ("coarsening"); }
      MPI_Barrier (MPI_COMM_WORLD);
            

      if (rank == 0) { tic (); }
      sprintf (filename, "palla_%5.5d", iframe++);
      tmsh.vtk_export (filename);
      if (rank == 0) { toc ("io"); }
      MPI_Barrier (MPI_COMM_WORLD);
  
    }
  
  if (rank == 0) {print_timing_report();}
  MPI_Barrier (MPI_COMM_WORLD);
  
  MPI_Finalize ();
  return 0;

}
