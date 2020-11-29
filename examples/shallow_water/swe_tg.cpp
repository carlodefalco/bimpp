/*
  Copyright (C) 2020 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <octave_file_io.h>

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

static int
refine_right_half (tmesh::quadrant_iterator quadrant)
{ return (quadrant->centroid (0) > 0.0 ? 1 : 0); }

// Re-Define tic and toc to add an MPI_Barrier
#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }

using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = std::vector<double>;         // Typedef for local q_0 vector

// Assemble vector from mesh.
// FIXME  the following two functions are copied over from
// "quad_operators.cpp" as they were not exported in an header,
// should find better way to avoid code duplication
void
assemble_vector (tmesh::quadrant_iterator& quadrant,
                 const std::array<double, 4>& locrhs,
                 Q1& rhs,
                 const ordering& ord = default_ord)
{

  std::vector<unsigned int> rows;
  rows.reserve (2);
  int i, r;

  for (i = 0; i < 4; ++i)
    {
      rows.clear ();

      if (! quadrant->is_hanging (i))
        rows.push_back (quadrant->gt (i));
      else
        {
          rows.push_back (quadrant->gparent (0, i));
          rows.push_back (quadrant->gparent (1, i));
        }

      for (r = 0; r < rows.size (); ++r)
        rhs[ord (rows[r])] +=
          locrhs[i] / rows.size ();
    }
}


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

  /*
  tmsh.set_refine_marker (refine_right_half);
  tmsh.refine (0);
  tmsh.set_refine_marker (refine_right_half);
  tmsh.refine (0);
  tmsh.set_refine_marker (refine_right_half);
  tmsh.refine (0);
  */
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

  Q1 Z    (ln_nodes);
  Z.get_owned_data ().assign (Z.get_owned_data ().size (), 0.0);
  
  peraire_stepper stp(tmsh, sol, ordh, ordUx, ordUy, Z);
  
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
           
            sol [ordh (quadrant->gt (ii))] = h0_fun (xx, yy);
            sol [ordUx(quadrant->gt (ii))] = Ux0_fun (xx, yy);
            sol [ordUy(quadrant->gt (ii))] = Uy0_fun (xx, yy);

            Z   [quadrant->gt (ii)] = Z_fun (xx, yy);
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

              Z   [quadrant->gparent(0,ii)] += 0.;
              Z   [quadrant->gparent(1,ii)] += 0.;

            }
        }
    }

  
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy);

  bim2a_solution_with_ghosts (tmsh, Z, replace_op);
  
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
  sprintf(filename, "swe_Z_0000");
  tmsh.octbin_export (filename, Z);

  std::vector<double> full_time_vector;
  full_time_vector.reserve (static_cast<int> (T/DELTAT));
  std::vector<double> save_time_vector;
  save_time_vector.reserve (static_cast<int> (T/SAVEDT));

  
  // Time loop
  double time = 0.0;
  double deltat = DELTAT;
  if(rank==0) {
    full_time_vector.push_back (0.0);
    save_time_vector.push_back (0.0);
  }
  int count = 0;
  
  double savecount = 0.0;
  
  while (time <= T)
    {      
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
          stp.update_halfstep ();
          stp.update_src ();
          stp.update_flux ();

          stp.update_state_incr ();
          assemble_vector (quadrant, stp.loc_incrh, incr, ordh);
          assemble_vector (quadrant, stp.loc_incrUx, incr, ordUx);
          assemble_vector (quadrant, stp.loc_incrUy, incr, ordUy);
        }

      incr.assemble ();
      TOC("Compute step");

      deltat = REDCDT * stp.get_dt ();
      if(rank==0)
        std::cout << "TIME = " << time << ", next dt = " << deltat << std::endl;
      
      MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&deltat), 1, MPI_INT, MPI_MIN, tmsh.comm);
      time += deltat;
      savecount += deltat;
      MPI_Bcast (static_cast<void*> (&time), 1, MPI_DOUBLE, 0, tmsh.comm);
      MPI_Bcast (static_cast<void*> (&savecount), 1, MPI_DOUBLE, 0, tmsh.comm);
      MPI_Barrier (tmsh.comm);

      // Print curent time
      if(rank==0)
        {
          std::cout << "TIME = " << time << ", last dt = " << deltat << std::endl;
          full_time_vector.push_back (time);
        }


      TIC();
      for (auto kk = 0; kk < incr.get_owned_data ().size (); kk++)
        sol.get_owned_data ()[kk] += deltat * incr.get_owned_data ()[kk] / mass.get_owned_data ()[kk];
      

      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep (); ++quadrant) {
        auto tree_idx = quadrant->get_tree_idx ();

        for (int i = 0; i < 4; ++i) {
          if (! quadrant->is_hanging (i)) {
              
            auto boundary_idx = quadrant->e (i);
            auto boundary_idxx = quadrant->ex (i);
            auto boundary_idxy = quadrant->ey (i);
            
            // If current node is on boundary 
            if (boundary_idx != tmesh::quadrant_t::NOT_ON_BOUNDARY) {

              // Loop over all the boundary conditions on h.
              for (size_t bc = 0; bc < bcsh.size (); ++bc)
                // If this boundary condition matches with
                // the current node.
                if (std::get<0> (bcsh[bc]) == tree_idx
                    && (std::get<1> (bcsh[bc]) == boundary_idx
                        || std::get<1> (bcsh[bc]) == boundary_idxx
                        || std::get<1> (bcsh[bc]) == boundary_idxy)) {                        
                  // Evaluate bc at current node.
                  sol[ordh (quadrant->gt (i))] = (std::get<2> (bcsh[bc])) (quadrant->p (0, i), quadrant->p (1, i));
                }

              // Loop over all the boundary conditions on Ux.
              for (size_t bc = 0; bc < bcsUx.size (); ++bc)
                // If this boundary condition matches with
                // the current node.
                if (std::get<0> (bcsUx[bc]) == tree_idx
                    && (std::get<1> (bcsUx[bc]) == boundary_idx
                        || std::get<1> (bcsUx[bc]) == boundary_idxx
                        || std::get<1> (bcsUx[bc]) == boundary_idxy)) {                    
                  // Evaluate bc at current node.
                  sol[ordUx (quadrant->gt (i))] = (std::get<2> (bcsUx[bc])) (quadrant->p (0, i), quadrant->p (1, i));
                }
              
              // Loop over all the boundary conditions on Uy.
              for (size_t bc = 0; bc < bcsUy.size (); ++bc)
                // If this boundary condition matches with
                // the current node.
                if (std::get<0> (bcsUy[bc]) == tree_idx
                    && (std::get<1> (bcsUy[bc]) == boundary_idx
                        || std::get<1> (bcsUy[bc]) == boundary_idxx
                        || std::get<1> (bcsUy[bc]) == boundary_idxy)) {                       
                  // Evaluate bc at current node.
                  sol[ordUy (quadrant->gt (i))] = (std::get<2> (bcsUy[bc])) (quadrant->p (0, i), quadrant->p (1, i));
                }
            }
          }
        }
      }
      sol.assemble (replace_op);            
      TOC("Apply increment");
      
      // Save solution      
      if (savecount >= SAVEDT) {
        TIC();
        if (rank == 0)
          std::cout << "savecount = " << savecount << std::endl;
        count++;
        save_time_vector.push_back (time);
        sprintf(filename, "swe_h_%4.4d",   count);
        tmsh.octbin_export (filename, sol, ordh);
        sprintf(filename, "swe_Ux_%4.4d",  count);
        tmsh.octbin_export (filename, sol, ordUx);
        sprintf(filename, "swe_Uy_%4.4d",  count);
        tmsh.octbin_export (filename, sol, ordUy);
        sprintf(filename, "swe_Z_%4.4d",  count);
        tmsh.octbin_export (filename, Z);
        savecount = 0.0;
        TOC("Exporting solution");
      }


    }

  if (rank == 0)
    {
      ColumnVector vtmp (save_time_vector.size ());
      std::copy (save_time_vector.begin (), save_time_vector.end (), vtmp.fortran_vec ());
      octave_io_mode m;
      octave_io_open ("timesteps.octbin", gz_write_mode, &m);
      octave_save ("save_time", vtmp);
      vtmp.resize (full_time_vector.size ());
      std::copy (full_time_vector.begin (), full_time_vector.end (), vtmp.fortran_vec ());
      octave_save ("full_time", vtmp);
      octave_io_close ();
    }

  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  MPI_Finalize ();

  return 0;
}
