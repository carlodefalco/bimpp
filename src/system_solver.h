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

constexpr int NUM_REFINEMENTS=5;

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; };

#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }




class system_solver
{
  private:
    int rank;
    int N;
    int NUM_ADAPT;
    int NUM_NON_ADAPT;
    int REFINE_ITER;
    double TOL_EST;
    double TOL_METRICS;
    double DELTAT;
    double T;
    double EPSU;
    double EPSV;
    double TAUU;
    double TAUV;
    double ALPHA;
    double BETA;
    double VBAR;
    double SIGMA;
    int count_progressive;
    int time_count;
    tmesh tmsh;
    std::vector<std::function<int(tmesh::idx_t)>> ord;
    std::vector<std::function<double(double x, double y)>> lambdas_init;
    std::vector<std::function<double(double x)>> lambdas_lin;
    std::vector<tmesh::idx_t> nnodes;
    std::vector<double> h_step;
    std::vector<double> estim;
    tmesh::idx_t gn_nodes;
    tmesh::idx_t ln_nodes;
    tmesh::idx_t ln_elements;
    int PRINT_PROG;
    int PRINT_SOL;
    int PRINT_GRAD;
    int PRINT_EST;





  public:
    void set_init_cond(distributed_vector* init);
    void solve();
    void print_sol_progressive( distributed_vector* to_print);
    void print_sol_analysis( distributed_vector* to_print, int m);
    void print_grad_analysis( std::vector<gradient<distributed_vector>>* gradients, int m);
    void print_estim_analysis( std::vector<double>* estim, int m);
    void obtain_global( distributed_vector* local, std::vector<distributed_vector>* global);
    void estimator_solution(std::vector<distributed_vector>* global, int m);
    void interpolation(distributed_vector* old_vec, distributed_vector* new_vec);
    void assemble_lin_sys (distributed_sparse_matrix* A, distributed_vector* rhs);
    void linear_solution(distributed_sparse_matrix* A, distributed_vector* rhs,
                    mumps* lin_solver, std::vector<double>* xa, std::vector<int>* ir, std::vector<int>* jc);
    void mesh_value_update();
    void mesh_init();
    void element_evaluate(std::vector<std::vector<double>>* cont_vec, std::vector<double>* coeff);
    void nodes_evaluate(std::vector<distributed_vector>* cont_vec,
                        std::vector<std::function<double(std::vector<double>)>>* lambdas_vec,
                        distributed_vector* uold);
    void assemble_matrix(distributed_sparse_matrix* A, distributed_vector* rhs, distributed_vector* uold);
    void solver_analysis(distributed_sparse_matrix* A, mumps * lin_solver, std::vector<double>* xa, std::vector<int> *ir, std::vector<int> *jc );
    void obtain_solution (mumps* lin_solver, distributed_vector* new_sol);

    system_solver(int rank_,
                  std::vector<std::function<double(double x, double y)>> lambdas_init_,
                  std::vector<std::function<double(double x)>> lambdas_lin_,
                  int N_=4,
                  int NUM_ADAPT_ = 20,
                  int NUM_NON_ADAPT_ = 4,
                  int REFINE_ITER_ = 2,
                  double TOL_EST_= 1e-4,
                  double TOL_METRICS_=1e5,
                  double DELTAT_ = 0.005,
                  double T_ = 5,
                  double EPSU_ = 0.05,
                  double EPSV_ = 0.05,
                  double TAUU_ = 1,
                  double TAUV_ = 100,
                  double ALPHA_ = 0.04,
                  double BETA_ = -0.9,
                  double VBAR_ = 0.0114559,
                  double SIGMA_ = 100,
                  int PRINT_PROG_=1,
                  int PRINT_SOL_=0,
                  int PRINT_GRAD_=0,
                  int PRINT_EST_=0):
      rank(rank_),
      N(N_),
      NUM_ADAPT(NUM_ADAPT_),
      NUM_NON_ADAPT(NUM_NON_ADAPT_),
      REFINE_ITER(REFINE_ITER_),
      TOL_EST(TOL_EST_),
      TOL_METRICS(TOL_METRICS_),
      DELTAT(DELTAT_),
      T(T_),
      EPSU(EPSU_),
      EPSV(EPSV_),
      TAUU(TAUU_),
      TAUV(TAUV_),
      ALPHA(ALPHA_),
      BETA(BETA_),
      VBAR(VBAR_),
      SIGMA(SIGMA_),
      lambdas_init(lambdas_init_),
      lambdas_lin(lambdas_lin_),
      time_count(0),
      count_progressive(0),
      PRINT_PROG(PRINT_PROG_),
      PRINT_SOL(PRINT_SOL_),
      PRINT_GRAD(PRINT_GRAD_),
      PRINT_EST(PRINT_EST_)
      {
        // Contruction of order function
        int n=N;
        auto lambda_gen = [n] (int i) {return [n,i](tmesh::idx_t gt){ return n*gt+i;};};
        for(int i = 0; i < N; i++)
        ord.push_back(lambda_gen(i));


      };

};
#endif /* SYSTEM_SOLVER_H */
