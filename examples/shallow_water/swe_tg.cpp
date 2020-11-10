/*
  Copyright (C) 2020 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>

#include <bim_distributed_vector.h>
#include <bim_timing.h>
#include <tmesh.h>
#include <quad_operators.h>

#include "swe_tg.h"
#include "stepper.h"

// Refinement rule
static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }


// Re-Define tic and toc to add an MPI_Barrier
#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }

int
main (int argc, char **argv)
{

  // Solution ordering
  ordering ordh  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 0> (gt); };
  ordering ordUx = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 1> (gt); };
  ordering ordUy = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 2> (gt); };

  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);


  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);

  tmesh::idx_t gn_nodes    = tmsh.num_global_nodes ();
  tmesh::idx_t ln_nodes    = tmsh.num_owned_nodes ();
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();

  Q1 sol  (ln_nodes * 3);
  Q1 incr (ln_nodes * 3);
  sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);
  incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);
  
  Q1 mass (ln_nodes * 3);
  bim2a_mass_vector (tmsh, mass, ordh);
  bim2a_mass_vector (tmsh, mass, ordUx);
  bim2a_mass_vector (tmsh, mass, ordUy);
  mass.assemble ();

  Q0 flux (ln_elements * 3);
  flux.assign (flux.size (), 0.0);
  
  peraire_stepper stp(tmsh, sol, ordh, ordUx, ordUy);
  
  // Buffer for export filename
  char filename[255]="";

  // Set initial value of state vector
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {      
      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){
            double xx=quadrant->p(0,ii);
            double yy=quadrant->p(1,ii);
           
            sol [ordh (quadrant->gt (ii))] = (8. - std::sin (pi * xx / 2. / 400.));
            sol [ordUx(quadrant->gt (ii))] = (8. - std::sin (pi * xx / 2. / 400.));
            sol [ordUy(quadrant->gt (ii))] = (8. - std::sin (pi * xx / 2. / 400.));
          }

          else
            {
              // touch parent nodes to set up distributed vector structure

              sol [ordh (quadrant->gparent(0,ii))] += 0.;
              sol [ordh (quadrant->gparent(1,ii))] += 0.;
              sol [ordUx(quadrant->gparent(0,ii))] += 0.;
              sol [ordUx(quadrant->gparent(1,ii))] += 0.;
              sol [ordUy(quadrant->gparent(0,ii))] += 0.;
              sol [ordUy(quadrant->gparent(1,ii))] += 0.;

            }
        }
    }

  
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy);
  
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUy);
  TOC ("compute initial condition");


  // Save initial conditions
  sprintf(filename, "swe_h_0000");
  tmsh.octbin_export (filename, sol, ordh);
  sprintf(filename, "swe_Ux_0000");
  tmsh.octbin_export (filename, sol, ordUx);
  sprintf(filename, "swe_Uy_0000");
  tmsh.octbin_export (filename, sol, ordUy);



  int count = 0;
  int savecount = 0;
  // Time loop
  for (double time = DELTAT; time <= T; time += DELTAT)
    {
      savecount++;
      
      // Print curent time
      if(rank==0)
        std::cout<<"TIME= "<<time<<std::endl;

      // Reset increment
      TIC();
      incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);
      incr.assemble (replace_op);
      TOC("Reset");

      TIC();
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
        {
          stp.set_quadrant (quadrant);
          stp.set_dt (DELTAT);
          stp.update_flux ();
          stp.get_flux (flux[ordh(quadrant->get_forest_quad_idx ())],
                        flux[ordUx(quadrant->get_forest_quad_idx ())],
                        flux[ordUy(quadrant->get_forest_quad_idx ())]);
        }
      incr.assemble ();
      TOC("Compute flux");
      
      TIC();
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
        {
          stp.set_quadrant (quadrant);
          stp.set_flux (flux[ordh(quadrant->get_forest_quad_idx ())],
                        flux[ordUx(quadrant->get_forest_quad_idx ())],
                        flux[ordUy(quadrant->get_forest_quad_idx ())]);
          stp.update_state ();
          assemble_vector (quadrant, stp.loc_incrh, incr, ordh);
          assemble_vector (quadrant, stp.loc_incrUx, incr, ordUx);
          assemble_vector (quadrant, stp.loc_incrUy, incr, ordUy);
        }
      incr.assemble ();
      TOC("Compute step");

      TIC();
      for (auto kk = 0; kk < incr.get_owned_data ().size (); kk++)
        sol.get_owned_data ()[kk] += DELTAT * incr.get_owned_data ()[kk] / mass.get_owned_data ()[kk];
      sol.assemble (replace_op);
      TOC("Apply increment");
      
      // Save solution
      TIC();
      if (savecount >= SKIPSAVE) {
        count++;
        sprintf(filename, "swe_h_%4.4d",   count);
        tmsh.octbin_export (filename, sol, ordh);
        sprintf(filename, "swe_Ux_%4.4d",  count);
        tmsh.octbin_export (filename, sol, ordUx);
        sprintf(filename, "swe_Uy_%4.4d",  count);
        tmsh.octbin_export (filename, sol, ordUy);
        savecount = 0;
      }
      TOC("Exporting solution");

    }


  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  MPI_Finalize ();

  return 0;
}
