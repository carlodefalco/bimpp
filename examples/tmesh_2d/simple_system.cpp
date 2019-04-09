/*
  Copyright (C) 2019 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <cmath>
#include <algorithm>

#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <tmesh.h>
#include <quad_operators.h>

constexpr int NUM_REFINEMENTS = 10;

constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {0., 0., 1., 0.,  1., 1., 0., 1.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {1, 2, 3, 4, 1};

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }

int
main (int argc, char **argv)
{

  using q1_vec  = q1_vec<distributed_vector>;

  ordering
    ord0 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<2, 0> (gt); },
    ord1 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<2, 1> (gt); };
    
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

  tmesh::idx_t gn_nodes = tmsh.num_global_nodes ();
  tmesh::idx_t ln_nodes = tmsh.num_owned_nodes ();
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();

  /// Allocate system matrix
  mumps *lin_solver = new mumps ();

  q1_vec u (ln_nodes * 2);
  u.get_owned_data ().assign (u.get_owned_data ().size (), 0.0);

  std::vector<double> xa;
  std::vector<int> ir, jc;
  
  distributed_sparse_matrix A;
  A.set_ranges (ln_nodes * 2);
  dirichlet_bcs bcs;
  
  for (int i = 0; i < 4; ++i)
    bcs.push_back (std::make_tuple (0, i, [](double x, double y){return .0;}));
      

  std::vector<double> ecoeff (ln_elements);
  q1_vec              ncoeff (ln_nodes);


  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      ecoeff[quadrant->get_forest_quad_idx ()] = 1.0;
      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii))
            ncoeff[quadrant->gt (ii)] = 1.0;
          else
            {
              ncoeff[quadrant->gparent (0, ii)] += 0.;
              ncoeff[quadrant->gparent (1, ii)] += 0.;
            }
        }
    }
  bim2a_solution_with_ghosts (tmsh, ncoeff, replace_op);
  TOC ("compute coefficient");


  TIC ();
  bim2a_advection_diffusion (tmsh, ecoeff, ncoeff, A, true, ord0, ord0);
  bim2a_advection_diffusion (tmsh, ecoeff, ncoeff, A, true, ord1, ord1);
  TOC ("assemble LHS");

  TIC ();
  bim2a_rhs (tmsh, ecoeff, ncoeff, u, ord0);
  bim2a_rhs (tmsh, ecoeff, ncoeff, u, ord1);
  TOC ("assemble RHS");
  
  
  TIC ();
  bim2a_dirichlet_bc (tmsh, bcs, A, u, ord0);
  bim2a_dirichlet_bc (tmsh, bcs, A, u, ord1);
  TOC ("apply BCS");
  
 
  TIC ();
  A.assemble ();
  u.assemble ();
  TOC ("communicate A and b");

  TIC ();
  lin_solver->set_lhs_distributed ();
  A.aij (xa, ir, jc, lin_solver->get_index_base ());
  lin_solver->set_distributed_lhs_structure (A.rows (), ir, jc);
  std::cout << "lin_solver->analyze () = "<< lin_solver->analyze () << std::endl;
  TOC ("solver analysis");


  TIC ();
  A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
  lin_solver->set_distributed_lhs_data (xa);
  TOC ("set LHS data");
  
  TIC ();
  lin_solver->set_rhs_distributed (u);
  TOC ("set RHS data");

  TIC ();  
  std::cout << "lin_solver->factorize () = " << lin_solver->factorize () << std::endl;
  TOC ("solver factorize");


  TIC ();  
  std::cout << "lin_solver->solve () = " << lin_solver->solve () << std::endl;
  TOC ("solver solve");

  // q1_vec result1(ln_nodes * 2);
  //{
    q1_vec result = lin_solver->get_distributed_solution ();
    // for (int ii = result1.get_range_start (); ii < result1.get_range_end (); ++ii)
    //   result1[ii] = result0[ii];
    bim2a_solution_with_ghosts (tmsh, result, replace_op, ord0, false);
    bim2a_solution_with_ghosts (tmsh, result, replace_op, ord1);
    // }
  
     tmsh.octbin_export (std::string ("simple_system_u_0").c_str (), result, ord0);
  tmsh.octbin_export (std::string ("simple_system_v_0").c_str (), result, ord1);

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  lin_solver->cleanup ();
  
  MPI_Finalize ();

  return 0;
}
