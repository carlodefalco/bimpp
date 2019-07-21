#ifndef SYSTEM_SETTINGS_H
#define SYSTEM_SETTINGS_H

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>
#include <functional>
#include <map>
#include "GetPot"
#include "muparser_fun.h"


class system_settings
{
  private:
    int N;                    // Number of equations
    int NUM_ADAPT;            // Number of total adaptation steps
    int NUM_NON_ADAPT;        // Number of steps done after an adaptation step
    int REFINE_ITER;          // Number of refinement steps in each adaptation step
    double DELTAT;            // Time step
    double TOL_EST;           // Toll for the estimator
    double TOL_METRICS;       // Toll for the metrics
    int PRINT_PROG;           // 1 to print progressively all the steps
    int PRINT_SOL;            // 1 to print solution in the form "solution_#idxsol_time_#refineiter_#processor".
                              // (refineiter != 0 only when a refinement step is done)
    int PRINT_GRAD;           // 1 to print the gradients
    int PRINT_EST;            // 1 to print the estimators



    // From file. Output as they are
    std::map<std::pair<int, int>, double> coeff_lap;                // Nonzero coefficients of the laplacian
    std::map<int, muparser_fun> fun_rhs;                            // Functions for Rhs
    std::map<int, muparser_fun> lambdas_init;                       // Initial conditions
    muparser_fun function_estimator;                                // Function of u and grad(u), to use as estimator.


    // Not from file. Need to be created as ouput
    std::map<std::pair<int, int>, muparser_fun> lambdas_vec_adv;    // Functions for the advection term (assembled)
    std::map<std::pair<int, int>, muparser_fun> lambdas_vec_reaz;   // Functions for the reaction term (assembled)
    std::map<int, muparser_fun> lambdas_vec_forc;                   // Functions for the forcing term (assembled)


    // Functions for reading from file
    void read_param(std::string const path,std::string const namevar,int const num_param,
      std::map<std::pair<int, int>, double> *output, GetPot fileData);
    void read_param(std::string const path,std::string const namevar,int const num_param,
      std::map<int, double> *output,GetPot fileData);
    void read_param(std::string const path,std::string const namevar,int const num_param,
      std::map<std::pair<int, int>, muparser_fun> *output, GetPot fileData);
    void read_param(std::string const path,std::string const namevar,int const num_param,
      std::map<int, muparser_fun> *output, GetPot fileData);


  public:

    friend class system_solver;

    // Read all the settings from a txt file "filename"
    void read_from_file_complete(std::string const filename);

    // Read only the equations parameters (coefficients and functions) from a txt file "filename"
    void read_from_file_reduced(std::string const filename);

    // Set the solver parameter
    void set_solver_parameter(int const N_,int const NUM_ADAPT_,int const NUM_NON_ADAPT_,int const REFINE_ITER_,
      double const DELTAT_,double const T_,int const PRINT_PROG_,int const PRINT_SOL_,int const PRINT_GRAD_,int const PRINT_EST_);

    // Set the tollerances for the mesh adaptivity
    void set_toll(double const TOL_EST_, double const TOL_METRICS_)
    {
      TOL_EST = TOL_EST_;
      TOL_METRICS = TOL_METRICS_;
    }

    // Default constructor
    system_settings():
      N(0),
      NUM_ADAPT(0),
      NUM_NON_ADAPT(0),
      REFINE_ITER(0),
      TOL_EST(0),
      TOL_METRICS(0),
      DELTAT(0),
      PRINT_PROG(0),
      PRINT_SOL(0),
      PRINT_GRAD(0),
      PRINT_EST(0)
      {};

    // Constructor for the solver parameters only
    system_settings(
                  int const N_,
                  int const NUM_ADAPT_,
                  int const NUM_NON_ADAPT_,
                  int const REFINE_ITER_,
                  double const TOL_EST_,
                  double const TOL_METRICS_,
                  double const DELTAT_,
                  double const T_,
                  int const PRINT_PROG_,
                  int const PRINT_SOL_,
                  int const PRINT_GRAD_,
                  int const PRINT_EST_):
      N(N_),
      NUM_ADAPT(NUM_ADAPT_),
      NUM_NON_ADAPT(NUM_NON_ADAPT_),
      REFINE_ITER(REFINE_ITER_),
      TOL_EST(TOL_EST_),
      TOL_METRICS(TOL_METRICS_),
      DELTAT(DELTAT_),
      PRINT_PROG(PRINT_PROG_),
      PRINT_SOL(PRINT_SOL_),
      PRINT_GRAD(PRINT_GRAD_),
      PRINT_EST(PRINT_EST_)
      {};

    // Constructor that reads all the settings from a txt file "filename"
    system_settings(std::string filename)
    {
      read_from_file_complete(filename);
    };

    // Constructor that reads the solver parameters and the equations parameters from a txt file "filename"
    system_settings(
                  int const N_,
                  int const NUM_ADAPT_,
                  int const NUM_NON_ADAPT_,
                  int const REFINE_ITER_,
                  double const TOL_EST_,
                  double const TOL_METRICS_,
                  double const DELTAT_,
                  double const T_,
                  int const PRINT_PROG_,
                  int const PRINT_SOL_,
                  int const PRINT_GRAD_,
                  int const PRINT_EST_,
                  std::string const filename):
    N(N_),
    NUM_ADAPT(NUM_ADAPT_),
    NUM_NON_ADAPT(NUM_NON_ADAPT_),
    REFINE_ITER(REFINE_ITER_),
    TOL_EST(TOL_EST_),
    TOL_METRICS(TOL_METRICS_),
    DELTAT(DELTAT_),
    PRINT_PROG(PRINT_PROG_),
    PRINT_SOL(PRINT_SOL_),
    PRINT_GRAD(PRINT_GRAD_),
    PRINT_EST(PRINT_EST_)
    {
      read_from_file_reduced(filename);
    };
};



#endif
