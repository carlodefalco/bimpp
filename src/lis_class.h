/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file lis_class.h
  \brief Interface for linear solver built with lis library.
*/

#ifndef HAVE_LIS_CLASS
#define HAVE_LIS_CLASS 1

#include <lis.h>
#include <string>
#include <bim_sparse.h>
#include <linear_solver.h>

/// Interface for linear solver built with lis library.
class lis: public linear_solver
{
private :

  LIS_SOLVER solver;
  LIS_MATRIX A;
  LIS_VECTOR b, x;
  LIS_INT iter;
  double time;

  LIS_INT *row,  *col;
  LIS_SCALAR *value;

  /// Stores rows of matrix in CSR format.
  std::vector<int> row_ptr;

  /// Stores columns of matrix in CSR format
  std::vector<int> jcol;

  /// Pointer to data of matrix.
  double *data;

  /// Pointer to values of rhs.
  double *rhs;

  /// Stores values of ordering map.
  std::vector<int> ordering_map;

  /// Pointer to values of initial guess.
  double *initial_guess;

  /// It's true if solver have initial guess passed by user.
  bool have_initial_guess;

  /// \brief Maximum number of iterations.
  /// \details [default = 1000].
  int max_iter;

  /// \brief Tolerance of linear solver.
  /// \details [default = 1e-12].
  double tolerance;

  /// \brief Name of iterative method.
  /// \details [default = bicg].
  std::string iterative_method;

  /// \brief Name of preconditioner.
  /// \details [default = none].
  std::string preconditioner;

  /// \brief Type of convergence condition.
  /// \details [default = norm 2 of residual]
  std::string convergence_condition;

  std::vector<int> map_i_s, map_row_s;
  std::vector<int> map_n, map_nnz;

  int i_s, row_s;
  int n, nnz, n_row;

  /// Index base used by specific linear solver.
  static const int index_base = 0;

  int rank, size;

public :

  /// \brief Option string.
  std::string option_string;

  /// Default costructor.
  lis (LIS_INT argc = 0, char * argv[] = NULL) :
    linear_solver ("LIS", "iterative")
  {
    lis_initialize (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    max_iter = 1000;
    tolerance = 1.0e-12;
    have_initial_guess = false;
    iterative_method = "bicg";
    preconditioner = "none";
    convergence_condition = "nrm2_r";

    option_string = "";
    
    lis_matrix_create (LIS_COMM_WORLD, &A);
    lis_vector_create (LIS_COMM_WORLD, &b);
    lis_vector_create (LIS_COMM_WORLD, &x);
    lis_solver_create (&solver);
  
  };

  /// Set-up the matrix structure.
  void
  set_lhs_structure
  (int n,
   std::vector<int> &ir,
   std::vector<int> &jc,
   matrix_format_t f = aij);

  /// Perform the analysis.
  int
  analyze ();

  /// Set matrix entries.
  void
  set_lhs_data (std::vector<double> &xa);

  /// Set the rhs.
  void
  set_rhs (std::vector<double> &rhs);

  /// Set the initial guess.
  void
  set_initial_guess (std::vector<double> &initial_guess);

  /// Prepare the solver.
  int
  factorize ();

  /// Solve the system.
  int
  solve ();

  /// Cleanup memory.
  void
  cleanup ();

  /// Set maximum number of iterations (default = 1000).
  void
  set_max_iterations (int max_iter_);

  /// Get maximum number of iterations.
  void
  get_max_iterations
  (int &max_iter_);

  /// Set tolerance of iterative method (default = 1.0e-12).
  void
  set_tolerance (double tol);

  /// Get tolerance of linear solver.
  void
  get_tolerance (double &tol);

  /// Set type of iterative method (default = cg).
  void
  set_iterative_method (const std::string &s);

  /// Get type of iterative method.
  void
  get_iterative_method (std::string &s);

  /// Set type of preconditioner (default = none).
  void
  set_preconditioner (const std::string &s);

  /// Get type of preconditioner.
  void
  get_preconditioner (std::string &s);

  /// Set convergence condition of iterative method
  /// (default = nrm2_r)
  void
  set_convergence_condition (const std::string &s);

  /// get convergence condition of iterative method
  void
  get_convergence_condition (std::string &s);

  /// LIS uses 0-based indexing
  inline int
  get_index_base ()
  { return index_base; }

};
#endif
