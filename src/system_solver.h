#ifndef SYSTEM_SOLVER_H
#define SYSTEM_SOLVER_H


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


#include "system_settings.h"

constexpr int NUM_REFINEMENTS=5;

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; };

#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }

class system_solver
{
  private:
    system_settings settings;                 // Problem settings

    int rank,size;

    int count_progressive;                    // Number of total
    int time_count;                           // Number of time interations

    tmesh tmsh;                               // Mesh

    std::vector<std::function<int(tmesh::idx_t)>> ord;    //  Vector of Neqations ordering

    std::vector<tmesh::idx_t> nnodes;
    std::vector<double> h_step;
    std::vector<double> estim;

    tmesh::idx_t gn_nodes;                    // Number of global nodes
    tmesh::idx_t ln_nodes;                    // Number of local nodes
    tmesh::idx_t ln_elements;                 // Number of local elements

    std::vector<int> index_grad_estim;        // Which gradients terms are present in the fuction used as estimator
    std::vector<int> index_sol_estim;         // Which solutions terms are present in the fuction used as estimator



    // Methods for the initializations
    void set_init_cond(distributed_vector* init);
    void mesh_value_update();
    void mesh_init();
    void reset_ordering(int const new_n)
    {
      ord.clear();
      int n=new_n;
      auto lambda_gen = [n] (int i) {return [n,i](tmesh::idx_t gt){ return n*gt+i;};};
      for(int i = 0; i < new_n; i++)
      ord.push_back(lambda_gen(i));
    }

    // Printing functions
    void print_sol_progressive( distributed_vector* to_print);
    void print_sol_analysis( distributed_vector* to_print, int m);
    void print_grad_analysis( std::vector<gradient<distributed_vector>>* gradients, int m);
    void print_estim_analysis( std::vector<double>* estim, int m);

    // Methods for assembling the linear system
    void element_evaluate(std::vector<double>* cont_vec, double value);
    void nodes_evaluate(distributed_vector* cont_vec,
                              muparser_fun* lambdas_vec,
                              distributed_vector* uold);
    void assemble_matrix(distributed_sparse_matrix* A, distributed_vector* rhs, distributed_vector* uold);

    // Methods for solving the linear system in the time loop
    void fake_assemble_matrix(distributed_sparse_matrix* A, distributed_vector* rhs);
    void assemble_lin_sys (distributed_sparse_matrix* A, distributed_vector* rhs);
    void linear_solution(distributed_sparse_matrix* A, distributed_vector* rhs,
                    mumps* lin_solver, std::vector<double>* xa, std::vector<int>* ir, std::vector<int>* jc);
    void solver_analysis(distributed_sparse_matrix* A, mumps * lin_solver,
                    std::vector<double>* xa, std::vector<int> *ir, std::vector<int> *jc );
    void obtain_solution (mumps* lin_solver, distributed_vector* new_sol);

    // Methods for the adaptivity in time
    void obtain_global( distributed_vector* local, std::vector<distributed_vector>* global);
    void estimator_solution(std::vector<distributed_vector>* global, int m);
    void interpolation(distributed_vector* old_vec, distributed_vector* new_vec);


  public:
    // Solve the problem loaded in settings
    void solve();

    // Set the problem from a system_settings object
    void set_problem(system_settings const new_problem)
    {
      settings = new_problem;
      reset_ordering(settings.N);
    }

    // Set the problem from a .txt file
    void set_problem(std::string const new_problem_file)
    {
      settings.read_from_file_complete(new_problem_file);
      reset_ordering(settings.N);
    }

    // Set the solver parameters and the equations parameters from a .txt file
    void set_problem(
                      int const N,
                      int const NUM_ADAPT,
                      int const NUM_NON_ADAPT,
                      int const REFINE_ITER,
                      double const TOL_EST,
                      double const TOL_METRICS,
                      double const DELTAT,
                      double const T,
                      int const PRINT_PROG,
                      int const PRINT_SOL,
                      int const PRINT_GRAD,
                      int const PRINT_EST,
                      std::string const new_problem_file)
    {
      settings.set_solver_parameter(N, NUM_ADAPT, NUM_NON_ADAPT, REFINE_ITER, DELTAT, T, PRINT_PROG, PRINT_SOL, PRINT_GRAD, PRINT_EST);
      settings.set_toll(TOL_EST,TOL_METRICS);

      settings.read_from_file_reduced(new_problem_file);


      reset_ordering(settings.N);
    }

    // Set the solver parameters
    void set_solver_parameter(
                          int N,
                          int NUM_ADAPT,
                          int NUM_NON_ADAPT,
                          int REFINE_ITER,
                          double TOL_EST,
                          double TOL_METRICS,
                          double DELTAT,
                          double T,
                          int PRINT_PROG,
                          int PRINT_SOL,
                          int PRINT_GRAD,
                          int PRINT_EST)
    {
      settings.set_solver_parameter(N, NUM_ADAPT, NUM_NON_ADAPT, REFINE_ITER, DELTAT, T, PRINT_PROG, PRINT_SOL, PRINT_GRAD, PRINT_EST);
      settings.set_toll(TOL_EST,TOL_METRICS);
      reset_ordering(settings.N);
    }

    // Default constructor
    system_solver():count_progressive(0),time_count(0),settings()
    {
      MPI_Comm_rank (MPI_COMM_WORLD, &rank);
      MPI_Comm_size (MPI_COMM_WORLD, &size);
    };

    // Constructor which reads the problem from a file .txt
    system_solver(std::string filename):count_progressive(0),time_count(0),settings(filename)
    {
      MPI_Comm_rank (MPI_COMM_WORLD, &rank);
      MPI_Comm_size (MPI_COMM_WORLD, &size);
      reset_ordering(settings.N);
    };

    // Constructor which reads the solver parameters and a file .txt containing the equations parameters
    system_solver(  int N_,
                    int NUM_ADAPT_,
                    int NUM_NON_ADAPT_,
                    int REFINE_ITER_,
                    double TOL_EST_,
                    double TOL_METRICS_,
                    double DELTAT_,
                    double T_,
                    int PRINT_PROG_,
                    int PRINT_SOL_,
                    int PRINT_GRAD_,
                    int PRINT_EST_,
                    std::string filename_):count_progressive(0),time_count(0),
                                            settings(N_,
                                                    NUM_ADAPT_,
                                                    NUM_NON_ADAPT_,
                                                    REFINE_ITER_,
                                                    TOL_EST_,
                                                    TOL_METRICS_,
                                                    DELTAT_,
                                                    T_,
                                                    PRINT_PROG_,
                                                    PRINT_SOL_,
                                                    PRINT_GRAD_,
                                                    PRINT_EST_,
                                                    filename_)
    {
      MPI_Comm_rank (MPI_COMM_WORLD, &rank);
      MPI_Comm_size (MPI_COMM_WORLD, &size);
      reset_ordering(settings.N);
    };


};
#endif
