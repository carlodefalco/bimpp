/*
  Copyright (C) 2019 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>



#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <bim_timing.h>
#include <mumps_class.h>
#include <tmesh.h>
#include <quad_operators.h>


// Setting parameters
constexpr int NUM_REFINEMENTS = 5;
constexpr int NUM_ADAPT = 40;
constexpr int NUM_NON_ADAPT = 10;
constexpr double EPS = 0.05;
constexpr double DELTAT = 0.005;
constexpr double T = 2;


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

  ////////////////INITIAL SETTING//////////////////////////////////////////////

  // Typedef for distributed vector
  using q1_vec  = q1_vec<distributed_vector>;


  // Manegement of solutions ordering
  ordering
    ord0 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<2, 0> (gt); },
    ord1 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<2, 1> (gt); };

  // Linearization techniques
  //OD2
  auto linc = [] (double x) -> double { return -1.5*x*x+0.5; };
  auto linf = [] (double x) -> double { return 0.5*x*(x*x+1); };

  //Eyre -> TO DO: add parameter tuning -beta || -x^3 + (beta+1)x
  //auto linc = [] (double x) -> double { return -2; };
  //auto linf = [] (double x) -> double { return (-x*(x*x - 3)); };

  //Linear splitting -> works with small time step
  // auto linc = [] (double x) -> double { return (1-x*x); };
  // auto linf = [] (double x) -> double { return 0; };

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

  // Allocate linear solver
  mumps *lin_solver = new mumps ();

  // Buffer for export filename
  char filename[255]="";


  // Squared epsilon parameter
  double eps2=EPS*EPS;


  int count = 0;

  std::vector<tmesh::idx_t> nnodes;
  std::vector<double> h_step;
  std::vector<double> estim;


///////////////////////////////FIRST TIME STEP/////////////////////////////////

  count++;

  // Print curent time
  if(rank==0)
  std::cout<<"TIME= "<<DELTAT*count<<std::endl;

  tmesh::idx_t gn_nodes = tmsh.num_global_nodes ();
  tmesh::idx_t ln_nodes = tmsh.num_owned_nodes ();
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();


  // Initial data
  q1_vec sold (ln_nodes * 2);
  sold.get_owned_data ().assign (sold.get_owned_data ().size (), 0.0);

  q1_vec sol (ln_nodes * 2);
  sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);


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
            sold[ord0(quadrant->gt (ii))] = std::sin(10*xx*yy);
            sol[ord0(quadrant->gt (ii))] = 0;
          }

          else
            {
              sold[ord0(quadrant->gparent(0,ii))] +=0.;
              sold[ord0(quadrant->gparent(1,ii))] +=0.;
              sol[ord0(quadrant->gparent(0,ii))] +=0.;
              sol[ord0(quadrant->gparent(1,ii))] +=0.;
            }
        }
    }
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord0, false);
  bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord1);
  TOC ("compute initial condition");

  // Save initial conditions
  sprintf(filename, "cahn_hilliard_u_0000");
  tmsh.octbin_export (filename, sold, ord0);


  // Containers construction
  std::vector<double> xa;
  std::vector<int> ir, jc;

  std::vector<double> lapcoeffu (ln_elements);
  std::vector<double> lapcoeffw (ln_elements);
  std::vector<double> reazuw (ln_elements);
  std::vector<double> reazwu (ln_elements);
  std::vector<double> ecoeff (ln_elements);

  q1_vec              ncoeff (ln_nodes);
  q1_vec              reazuu (ln_nodes);
  q1_vec              fu (ln_nodes);
  q1_vec              fw (ln_nodes);

  distributed_sparse_matrix A;
  A.set_ranges (ln_nodes * 2);


  // Initialize parameters
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      lapcoeffu[quadrant->get_forest_quad_idx ()] = eps2;
      lapcoeffw[quadrant->get_forest_quad_idx ()] = DELTAT;
      reazuw[quadrant->get_forest_quad_idx ()] = 1.0;
      reazwu[quadrant->get_forest_quad_idx ()] = -1.0;
      ecoeff[quadrant->get_forest_quad_idx ()] = 1.0;

      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){
            ncoeff[quadrant->gt (ii)] = 1.0;
            reazuu[quadrant->gt (ii)] = -linc(sold[ord0(quadrant->gt (ii))]);
            fu[quadrant->gt (ii)] = linf(sold[ord0(quadrant->gt (ii))]);
            fw[quadrant->gt (ii)] = -sold[ord0(quadrant->gt (ii))];
          }

          else
            {
              ncoeff[quadrant->gparent (0, ii)] += 0.;
              ncoeff[quadrant->gparent (1, ii)] += 0.;
              reazuu[quadrant->gparent (0, ii)] += 0.;
              reazuu[quadrant->gparent (1, ii)] += 0.;
              fu[quadrant->gparent (0, ii)] += 0.;
              fu[quadrant->gparent (1, ii)] += 0.;
              fw[quadrant->gparent (0, ii)] += 0.;
              fw[quadrant->gparent (1, ii)] += 0.;
            }
        }
    }
  ncoeff.assemble (replace_op);
  reazuu.assemble (replace_op);
  fu.assemble (replace_op);
  fw.assemble (replace_op);
  TOC ("compute coefficients");


  // Matrix construction
  TIC ();
  bim2a_laplacian(tmsh, lapcoeffu, A, ord0, ord0);
  bim2a_laplacian(tmsh, lapcoeffw, A, ord1, ord1);

  bim2a_reaction(tmsh, ecoeff, reazuu, A, ord0, ord0);
  bim2a_reaction(tmsh, reazuw, ncoeff, A, ord0, ord1);
  bim2a_reaction(tmsh, reazwu, ncoeff, A, ord1, ord0);
  A.assemble();
  TOC ("assemble LHS");


  // RHS construction
  TIC();
  bim2a_rhs (tmsh, ecoeff, fu, sol, ord0);
  bim2a_rhs (tmsh, ecoeff, fw, sol, ord1);
  sol.assemble();
  TOC ("assemble RHS");


  // Solver analysis
  TIC ();
  lin_solver->set_lhs_distributed ();
  A.aij (xa, ir, jc, lin_solver->get_index_base ());
  lin_solver->set_distributed_lhs_structure (A.rows (), ir, jc);
  std::cout << "lin_solver->analyze () = "<< lin_solver->analyze () << std::endl;
  TOC ("solver analysis");


  // Matrix update
  TIC ();
  A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
  lin_solver->set_distributed_lhs_data (xa);
  TOC ("set LHS data");


  // Factorization
  TIC ();
  std::cout << "lin_solver->factorize () = " << lin_solver->factorize () << std::endl;
  TOC ("solver factorize");


  // Set RHS data
  TIC ();
  lin_solver->set_rhs_distributed (sol);
  TOC ("set RHS data");


  // Solution
  TIC ();
  std::cout << "lin_solver->solve () = " << lin_solver->solve () << std::endl;
  TOC ("solver solve");


  // Copy solution
  TIC();
  q1_vec result = lin_solver->get_distributed_solution ();
  for (int idx = sold.get_range_start (); idx < sold.get_range_end (); ++idx)
    sold (idx) = result (idx);
  sold.assemble (replace_op);
  TOC("Obtaining solution");


  // Save solution
  TIC();
  sprintf(filename, "cahn_hilliard_u_%4.4d",count);
  tmsh.octbin_export (filename, sold, ord0);
  TOC("Exporting solution");

  lin_solver->cleanup ();


  /*
  PROBLEMI:

  Sono commentate le righe
  519
  520
  521
  598-602

  Sono le uniche righe dove, dopo il raffinamento, si accede a sold.
  Il problema è quindi che sold non "parla" con la nuova mesh.
  Ho provato a ridichiararla da 0 sulla nuova mesh ma non compila.
  Se invece si toglie il raffinamento va tutto.
  */



////////////////////////////////////////TIME CICLE//////////////////////////////

  for( int i=0 ; i< NUM_ADAPT ; i++){

    //////////////////////////////ESTIMATOR////////////////////////////////////////


        TIC();
        q1_vec only_u(gn_nodes);
        for (int idx = 0; idx < gn_nodes; ++idx)
          only_u(idx)=result(ord0(idx));
        only_u.assemble(replace_op);
        TOC("copy only u");



        TIC();
         double tol = 1e-4;
         gradient<q1_vec> du = bim2c_quadtree_pde_recovered_gradient(tmsh, only_u);
         q2_vec u_star = bim2c_quadtree_pde_recovered_solution(tmsh, only_u, du);
        TOC("gradient and ustar");

         sprintf(filename, "du_x_%4.4d",i+1);
         tmsh.octbin_export (filename, du.first);

         sprintf(filename, "du_y_%4.4d",i+1);
         tmsh.octbin_export (filename, du.second);


         auto estimator = [& u_star, & only_u] (tmesh::quadrant_iterator q)
           { return estimator_sol (q, u_star, only_u); };


         tmsh.set_metrics_marker (estimator, tol, 4, 2, 2);


         // Compute metrics and h.
         std::vector<double> metrics(ln_elements);

         double hx = 0, hy = 0,
         h = std::numeric_limits<double>::max (),
         global_h = 0;
         double est = 0, global_est = 0;

         for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
               ++quadrant)
             {
               metrics[quadrant->get_forest_quad_idx ()] =
             estimator(quadrant)* std::sqrt (tmsh.num_global_quadrants ())
               / tol;

               hx = quadrant->p(0, 1) - quadrant->p(0, 0);
               hy = quadrant->p(1, 2) - quadrant->p(1, 0);

               h = std::min(h, std::sqrt(hx*hx + hy*hy));

               est += std::pow(estimator(quadrant), 2);
             }

             sprintf(filename, "metrics_hx_%4.4d",i+1);
             tmsh.octbin_export_quadrant (filename, metrics);

             MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
             MPI_Reduce(&est, &global_est, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
             global_est = std::sqrt(global_est);

             nnodes.push_back (tmsh.num_global_nodes ());
             h_step.push_back (global_h);
             estim.push_back (global_est);

           TIC();
              tmsh.metrics_refine (1e5);  // RAFFINAMENTO
          TOC("refine");


    //////////////////////////////////////////////////////////////////////////////

            // MI PREPARO A FARE 10 PASSI TEMPORALI

            mumps *lin_solver = new mumps ();

            // Ottengo i parametri della mesh corrente
            TIC();
            gn_nodes = tmsh.num_global_nodes ();
            ln_nodes = tmsh.num_owned_nodes ();
            ln_elements = tmsh.num_local_quadrants ();
            TOC("Obtaining new parameters");


            // Interpolo sold sulla nuova mesh
            TIC();
            q1_vec new_result(gn_nodes*2);
            interpolate_vector (tmsh, result, new_result, ord0);
            q1_vec sold(ln_nodes*2);
            sold.get_owned_data ().assign (sold.get_owned_data ().size (), 0.0);
            for (int idx = sold.get_range_start (); idx < sold.get_range_end (); ++idx)
              sold (idx) = new_result (idx);
            sold.assemble (replace_op);
            TOC("Interpolation");


            TIC();
            distributed_sparse_matrix A;
            A.set_ranges (ln_nodes * 2);

            q1_vec sol (ln_nodes * 2);
            sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);

            std::vector<double> xa;
            std::vector<int> ir, jc;

            std::vector<double> lapcoeffu (ln_elements);
            std::vector<double> lapcoeffw (ln_elements);
            std::vector<double> reazuw (ln_elements);
            std::vector<double> reazwu (ln_elements);
            std::vector<double> ecoeff (ln_elements);

            q1_vec              ncoeff (ln_nodes);
            q1_vec              reazuu (ln_nodes);
            q1_vec              fu (ln_nodes);
            q1_vec              fw (ln_nodes);
            TOC("Containers recontruction")

            // Compute constant (in time ) coeffficients
            TIC ();
            for (auto quadrant = tmsh.begin_quadrant_sweep ();
                 quadrant != tmsh.end_quadrant_sweep ();
                 ++quadrant)
              {
                lapcoeffu[quadrant->get_forest_quad_idx ()] = eps2;
                lapcoeffw[quadrant->get_forest_quad_idx ()] = DELTAT;
                reazuw[quadrant->get_forest_quad_idx ()] = 1.0;
                reazwu[quadrant->get_forest_quad_idx ()] = -1.0;
                ecoeff[quadrant->get_forest_quad_idx ()] = 1.0;

                for (int ii = 0; ii < 4; ++ii)
                  {
                    if (! quadrant->is_hanging (ii)){
                      ncoeff[quadrant->gt (ii)] = 1.0;
                    }

                    else
                      {
                        ncoeff[quadrant->gparent (0, ii)] += 0.;
                        ncoeff[quadrant->gparent (1, ii)] += 0.;
                      }
                  }
              }
            ncoeff.assemble (replace_op);
            TOC ("compute coefficients");


            // Matrix and RHS construction (fake, just to set solver)
            TIC ();
            bim2a_laplacian(tmsh, ecoeff, A, ord0, ord0);
            bim2a_laplacian(tmsh, ecoeff, A, ord1, ord1);

            bim2a_reaction(tmsh, ecoeff, ncoeff, A, ord0, ord1);
            bim2a_reaction(tmsh, ecoeff, ncoeff, A, ord1, ord0);
            A.assemble ();
            TOC ("assemble LHS");

            TIC();
            bim2a_rhs (tmsh, ecoeff, ncoeff, sol, ord0);
            bim2a_rhs (tmsh, ecoeff, ncoeff, sol, ord1);
            sol.assemble ();
            TOC ("assemble RHS");


            // Solver analysis
            TIC ();
            lin_solver->set_lhs_distributed ();
            A.aij (xa, ir, jc, lin_solver->get_index_base ());
            lin_solver->set_distributed_lhs_structure (A.rows (), ir, jc);
            std::cout << "lin_solver->analyze () = "<< lin_solver->analyze () << std::endl;
            TOC ("solver analysis");



            for (int j =0 ; j< NUM_NON_ADAPT; j++){
              count++;

              // Print curent time
              if(rank==0)
              std::cout<<"TIME= "<<DELTAT*count<<std::endl;


             // Reset containers
             TIC();
             A.reset ();
             sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);
             sol.assemble (replace_op);
             TOC("Resetting")


              // Initialize non constant (in time) parameters
              TIC();
              for (auto quadrant = tmsh.begin_quadrant_sweep ();
                 quadrant != tmsh.end_quadrant_sweep ();
                 ++quadrant){
                   for (int ii = 0; ii < 4; ++ii)
                   {
                     if (! quadrant->is_hanging (ii)){
                       ncoeff[quadrant->gt (ii)] = 1.0;
                       reazuu[quadrant->gt (ii)] = -linc(sold[ord0(quadrant->gt (ii))]);
                       fu[quadrant->gt (ii)] = linf(sold[ord0(quadrant->gt (ii))]);
                       fw[quadrant->gt (ii)] = -sold[ord0(quadrant->gt (ii))];
                     }
                     else
                     {
                       reazuu[quadrant->gparent (0, ii)] += 0.;
                       reazuu[quadrant->gparent (1, ii)] += 0.;
                       fu[quadrant->gparent (0, ii)] += 0.;
                       fu[quadrant->gparent (1, ii)] += 0.;
                       fw[quadrant->gparent (0, ii)] += 0.;
                       fw[quadrant->gparent (1, ii)] += 0.;
                     }
                   }
              }
              reazuu.assemble (replace_op);
              fu.assemble (replace_op);
              fw.assemble (replace_op);
              TOC("Update coefficients");


              // Matrix construction
              TIC ();
              bim2a_laplacian(tmsh, lapcoeffu, A, ord0, ord0);
              bim2a_laplacian(tmsh, lapcoeffw, A, ord1, ord1);

              bim2a_reaction(tmsh, ecoeff, reazuu, A, ord0, ord0);
              bim2a_reaction(tmsh, reazuw, ncoeff, A, ord0, ord1);
              bim2a_reaction(tmsh, reazwu, ncoeff, A, ord1, ord0);
              TOC ("assemble LHS");

              // RHS construction
              TIC();
              bim2a_rhs (tmsh, ecoeff, fu, sol, ord0);
              bim2a_rhs (tmsh, ecoeff, fw, sol, ord1);
              TOC ("assemble RHS");


              // Communicate matrix and RHS
              TIC ();
              A.assemble ();
              sol.assemble ();
              TOC ("communicate A and b");


              // Matrix update
              TIC ();
              A.aij_update (xa, ir, jc, lin_solver->get_index_base ());
              lin_solver->set_distributed_lhs_data (xa);
              TOC ("set LHS data");


              // Factorization
              TIC ();
              std::cout << "lin_solver->factorize () = " << lin_solver->factorize () << std::endl;
              TOC ("solver factorize");


              // Set RHS data
              TIC ();
              lin_solver->set_rhs_distributed (sol);
              TOC ("set RHS data");


              // Solution
              TIC ();
              std::cout << "lin_solver->solve () = " << lin_solver->solve () << std::endl;
              TOC ("solver solve");


              // Copy solution
              TIC();
              result = lin_solver->get_distributed_solution ();
              for (int idx = sold.get_range_start (); idx < sold.get_range_end (); ++idx)
                sold (idx) = result (idx);
              bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord0, false);
              bim2a_solution_with_ghosts (tmsh, sold, replace_op, ord1);
              TOC("Obtaining solution");


              TIC();
              sprintf(filename, "cahn_hilliard_u_%4.4d",count);
              tmsh.octbin_export (filename, sold, ord0);
              TOC("Exporting solution");
            }
            lin_solver->cleanup ();

    }




    // Print number of nodes and hx
    if (rank == 0)
      for (unsigned step = 0; step < nnodes.size(); ++step)
        std::cout << "Step " << step << ", #nodes: "
                  << nnodes[step] << ", h: "
                  << h_step[step] << ", estimator: "
  		<< estim[step] << std::endl;


  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  MPI_Finalize ();


  return 0;
}
