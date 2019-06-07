/*
  Copyright (C) 2019 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include <time.h>


#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <tmesh.h>
#include <quad_operators.h>

// Setting parameters
constexpr int NUM_REFINEMENTS = 4;
constexpr double DELTAT = 0.005;
constexpr double T = 10;
constexpr double EPSU = 0.05;
constexpr double EPSV = 0.05;
constexpr double TAUU = 1;
constexpr double TAUV = 100;
constexpr double ALPHA = 0.04;
constexpr double BETA = -0.9;
constexpr double VBAR = 0.0114559;
constexpr double SIGMA = 100;



// Connectivity of local element
constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {0., 0., 1., 0.,  1., 1., 0., 1.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {1, 2, 3, 4, 1};

// Refinement rule
static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

// Define tic and toc
#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }

int
main (int argc, char **argv)
{

  using q1_vec  = q1_vec<distributed_vector>;                                       // Typedef for distributed vector

// Manegement of solutions ordering
  ordering
    ord0 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<4, 0> (gt); },
    ord1 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<4, 1> (gt); },
    ord2 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<4, 2> (gt); },
    ord3 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<4, 3> (gt); };

// Linearization techniques
  auto linc = [] (double x) -> double { return -1.5*x*x+0.5; };
  auto linf = [] (double x) -> double { return 0.5*x*(x*x+1); };



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

  tmesh::idx_t gn_nodes = tmsh.num_global_nodes ();
  tmesh::idx_t ln_nodes = tmsh.num_owned_nodes ();
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();



  // Allocate linear solver
  mumps *lin_solver = new mumps ();

  // Allocate initial data container
  q1_vec sold (ln_nodes * 4);
  sold.get_owned_data ().assign (sold.get_owned_data ().size (), 0.0);


  q1_vec sol (ln_nodes * 4);


  // Allocation of vector to use sparse matrix
  std::vector<double> xa;
  std::vector<int> ir, jc;

  // Allocation of container for system coefficients
  std::vector<double> lapcoeffu (ln_elements);
  q1_vec              reazuu (ln_nodes);
  std::vector<double> ones (ln_elements);
  q1_vec              reazuv (ln_nodes);
  q1_vec              fu (ln_nodes);

  std::vector<double> times (ln_elements);
  std::vector<double> reazwuu (ln_elements);
  q1_vec              fwu (ln_nodes);

  std::vector<double> lapcoeffv (ln_elements);
  std::vector<double> reazvu (ln_elements);
  q1_vec              reazvv (ln_nodes);
  q1_vec              fv (ln_nodes);

  std::vector<double> reazwvv (ln_elements);
  q1_vec              fwv (ln_nodes);

  q1_vec              ncoeff (ln_nodes);

  // Buffer for export filename
  char filename[255]="";

  // Squared epsilon parameters
  double epsu2=EPSU*EPSU;
  double epsv2=EPSV*EPSV;

  // Initialize random number generator
  std::default_random_engine generator(1+rank*10*clock());
  std::uniform_real_distribution<double> distribution(-0.5,0.5);




  // Initialize constant (in time) parameters and initial data
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      lapcoeffu[quadrant->get_forest_quad_idx ()] = epsu2;
      ones[quadrant->get_forest_quad_idx ()] = 1.0;
      times[quadrant->get_forest_quad_idx ()] = DELTAT;
      reazwuu[quadrant->get_forest_quad_idx ()] = -TAUU;
      lapcoeffv[quadrant->get_forest_quad_idx ()] = epsv2;
      reazvu[quadrant->get_forest_quad_idx ()] = ALPHA;
      reazwvv[quadrant->get_forest_quad_idx ()] = -(SIGMA*DELTAT+TAUV);

      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){
            ncoeff[quadrant->gt (ii)] = 1.0;
            double xx=quadrant->p(0,ii);
            double yy=quadrant->p(1,ii);
            sold[ord0(quadrant->gt (ii))] = std::sin(10*xx*yy);//distribution(generator);
            sold[ord2(quadrant->gt (ii))] = std::cos(10*(xx-yy))*xx*yy; //distribution(generator);
          }

          else
            {
              ncoeff[quadrant->gparent (0, ii)] += 0.;
              ncoeff[quadrant->gparent (1, ii)] += 0.;
              sold[ord0(quadrant->gparent(0,ii))] +=0.;
              sold[ord0(quadrant->gparent(1,ii))] +=0.;
              sold[ord2(quadrant->gparent(0,ii))] +=0.;
              sold[ord2(quadrant->gparent(1,ii))] +=0.;
            }
        }
    }

  bim2a_solution_with_ghosts (tmsh, ncoeff, replace_op);
  TOC ("compute coefficient");


  // Save initial conditions
  sprintf(filename, "cahn_hilliard_u_0000");
  tmsh.octbin_export (filename, sold, ord0);
  sprintf(filename, "cahn_hilliard_v_0000");
  tmsh.octbin_export (filename, sold, ord2);

  distributed_sparse_matrix A;
  A.set_ranges (ln_nodes * 4);
  
  // Time cycle
  for( int count=1; count <= T/DELTAT; count++){

    // Print current time
    if(rank==0)
    std::cout<<"TIME= "<<count*DELTAT<<std::endl;



    // Reset containers for linear solver -> must be improved
    TIC();   
    A.reset ();
    
    sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);
    sol.assemble (replace_op));
  
    TOC("Resetting");


    // Initialize non constant (in time) parameters
        TIC ();
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant){
               for (int ii = 0; ii < 4; ++ii)
               {
                 if (! quadrant->is_hanging (ii)){
                   reazuu[quadrant->gt (ii)] = -linc(sold[ord0(quadrant->gt (ii))]);
                   reazuv[quadrant->gt (ii)] = ALPHA+BETA*sold[ord2(quadrant->gt (ii))];
                   fu[quadrant->gt (ii)] = linf(sold[ord0(quadrant->gt (ii))]);
                   fwu[quadrant->gt (ii)] = -TAUU*sold[ord0(quadrant->gt (ii))];
                   reazvv[quadrant->gt (ii)] = 2*BETA*sold[ord0(quadrant->gt (ii))]-linc(sold[ord2(quadrant->gt (ii))]);
                   fv[quadrant->gt (ii)] = linf(sold[ord2(quadrant->gt (ii))]);
                   fwv[quadrant->gt (ii)] = -SIGMA*DELTAT*VBAR-TAUV*sold[ord2(quadrant->gt (ii))];
                 }
                 else
                 {
                   reazuu[quadrant->gparent (0, ii)] += 0.;
                   reazuu[quadrant->gparent (1, ii)] += 0.;
                   reazuv[quadrant->gparent (0, ii)] += 0.;
                   reazuv[quadrant->gparent (1, ii)] += 0.;
                   fu[quadrant->gparent (0, ii)] += 0.;
                   fu[quadrant->gparent (1, ii)] += 0.;
                   fwu[quadrant->gparent (0, ii)] += 0.;
                   fwu[quadrant->gparent (1, ii)] += 0.;
                   reazvv[quadrant->gparent (0, ii)] += 0.;
                   reazvv[quadrant->gparent (1, ii)] += 0.;
                   fv[quadrant->gparent (0, ii)] += 0.;
                   fv[quadrant->gparent (1, ii)] += 0.;
                   fwv[quadrant->gparent (0, ii)] += 0.;
                   fwv[quadrant->gparent (1, ii)] += 0.;
                 }
               }
          }

          // Is it necessary ?
          bim2a_solution_with_ghosts (tmsh, reazuu, replace_op);
          bim2a_solution_with_ghosts (tmsh, reazuv, replace_op);
          bim2a_solution_with_ghosts (tmsh, fu, replace_op);
          bim2a_solution_with_ghosts (tmsh, fwu, replace_op);
          bim2a_solution_with_ghosts (tmsh, reazvv, replace_op);
          bim2a_solution_with_ghosts (tmsh, fv, replace_op);
          bim2a_solution_with_ghosts (tmsh, fwv, replace_op);
          TOC("Update coefficients");

          // Matrix construction
          TIC ();
          bim2a_laplacian(tmsh, lapcoeffu, A, ord0, ord0);
          bim2a_reaction(tmsh, ones, reazuu, A, ord0, ord0);
          bim2a_reaction(tmsh, ones, ncoeff, A, ord0, ord1); //reazuwu
          bim2a_reaction(tmsh, ones, reazuv, A, ord0, ord2);


          bim2a_laplacian(tmsh, times, A, ord1, ord1);            //lapcoeffwu
          bim2a_reaction(tmsh, reazwuu, ncoeff, A, ord1, ord0);


          bim2a_laplacian(tmsh, lapcoeffv, A, ord2, ord2);
          bim2a_reaction(tmsh, reazvu, ncoeff, A, ord2, ord0);
          bim2a_reaction(tmsh, ones, reazvv, A, ord2, ord2);
          bim2a_reaction(tmsh, ones, ncoeff, A, ord2, ord3); //reazvwv


          bim2a_laplacian(tmsh, times, A, ord3, ord3);          //lapcoeffwv
          bim2a_reaction(tmsh, reazwvv, ncoeff, A, ord3, ord2);
          TOC ("assemble LHS");



          TIC();
          bim2a_rhs (tmsh, ones, fu, sol, ord0);
          bim2a_rhs (tmsh, ones, fwu, sol, ord1);
          bim2a_rhs (tmsh, ones, fv, sol, ord2);
          bim2a_rhs (tmsh, ones, fwv, sol, ord3);
          TOC ("assemble RHS");


          // Communicate matrix and RHS
          TIC ();
          A.assemble ();
          sol.assemble ();
          TOC ("communicate A and b");


          // Solver analysis
          TIC ();
          lin_solver->set_lhs_distributed ();
          A.aij (xa, ir, jc, lin_solver->get_index_base ());
          lin_solver->set_distributed_lhs_structure (A.rows (), ir, jc);
          std::cout << "lin_solver->analyze () = "<< lin_solver->analyze () << std::endl;
          TOC ("solver analysis");


          // Set LHS data
          TIC ();
          A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
          lin_solver->set_distributed_lhs_data (xa);
          TOC ("set LHS data");


          // Set RHS data
          TIC ();
          lin_solver->set_rhs_distributed (sol);
          TOC ("set RHS data");


          // Factorization
          TIC ();
          std::cout << "lin_solver->factorize () = " << lin_solver->factorize () << std::endl;
          TOC ("solver factorize");


          // Solution
          TIC ();
          std::cout << "lin_solver->solve () = " << lin_solver->solve () << std::endl;
          TOC ("solver solve");


          // Copy solution
          TIC();
          q1_vec result = lin_solver->get_distributed_solution ();
          bim2a_solution_with_ghosts (tmsh, result, replace_op, ord0, false);
          bim2a_solution_with_ghosts (tmsh, result, replace_op, ord2);
          bim2a_solution_with_ghosts (tmsh, result, replace_op, ord1);
          bim2a_solution_with_ghosts (tmsh, result, replace_op, ord3);
          TOC("Obtaining solution");


          // Save solution
          TIC();
          sprintf(filename, "cahn_hilliard_u_%4.4d",count);
          tmsh.octbin_export (filename, result, ord0);
          sprintf(filename, "cahn_hilliard_v_%4.4d",count);
          tmsh.octbin_export (filename, result, ord2);
          TOC("Exporting solution");


          // Update old solution
          TIC();
          sold=result;
          TOC("Updating solution");
}

  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  // Clean linear solver
  lin_solver->cleanup ();

  MPI_Finalize ();

  return 0;
}
