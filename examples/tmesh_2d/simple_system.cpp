/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <cmath>
#include <algorithm>

#include <bim_timing.h>
#include <mumps_class.h>
#include <tmesh.h>
#include <quad_operators.h>


constexpr int NUM_REFINEMENTS = 4;

constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {0., 0., 1., 0.,  1., 1., 0., 1.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {1, 2, 3, 4, 1};

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

int
main (int argc, char **argv)
{

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

  tmesh::idx_t n_nodes = tmsh.num_global_nodes ();
  tmesh::idx_t n_elements = tmsh.num_local_quadrants ();

  /// Allocate system matrix
  mumps *lin_solver = new mumps ();

  std::vector<double> u;
  u.assign (n_nodes*2, 0.0);

  std::vector<double> xa;
  std::vector<int> ir, jc;
  
  sparse_matrix A;
  A.resize (n_nodes*2);
  dirichlet_bcs bcs;
  
  for (int i = 0; i < 4; ++i)
    bcs.push_back (std::make_tuple (0, i, [](double x, double y){return .0;}));
      
  bim2a_structure (tmsh, A, ord0);
  bim2a_structure (tmsh, A, ord1);
  
  A.aij (xa, ir, jc, lin_solver->get_index_base ());

  std::vector<double> ecoeff (n_nodes);
  std::vector<double> ncoeff (n_elements);

  ecoeff.assign (n_elements, 1.0);
  ncoeff.assign (n_nodes, 1.0);

  lin_solver->set_lhs_distributed ();
  lin_solver->set_distributed_lhs_structure (A.rows (), ir, jc);
  lin_solver->analyze ();

  {
    if (rank == 0) tic ();
    A.reset ();
             
    bim2a_advection_diffusion (tmsh, ecoeff, ncoeff, A, ord0, ord0);
    bim2a_advection_diffusion (tmsh, ecoeff, ncoeff, A, ord1, ord1);
    bim2a_rhs (tmsh, ecoeff, ncoeff, u, ord0, ord0);
    bim2a_rhs (tmsh, ecoeff, ncoeff, u, ord1, ord0);

    
    bim2a_dirichlet_bc (tmsh, bcs, A, u, ord0);
    bim2a_dirichlet_bc (tmsh, bcs, A, u, ord1);
    
    std::cout << A << std::endl;
 
    MPI_Barrier (MPI_COMM_WORLD);
    if (rank == 0) toc ("assembly");
  }
  
  {
    if (rank == 0) tic ();
    A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
    lin_solver->set_distributed_lhs_data (xa);
    lin_solver->set_rhs (u);
    
    lin_solver->factorize ();
    lin_solver->solve ();
    MPI_Barrier (MPI_COMM_WORLD);
    if (rank == 0) toc ("solve");
  }
  
  tmsh.octbin_export (std::string ("simple_u").c_str (), u);

  print_timing_report ();
  lin_solver->cleanup ();
  MPI_Finalize ();

  return 0;
}



