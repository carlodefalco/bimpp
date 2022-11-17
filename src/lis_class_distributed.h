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

#include <iostream>
#include <lis.h>
#include <string>
#include <vector>
#include <linear_solver.h>


/// Interface for linear solver built with lis library.
class lis_distributed: public linear_solver
{
private :

  LIS_SOLVER solver;
  LIS_MATRIX A;
  LIS_VECTOR b, x;
  LIS_INT iter;
  double time;
  bool initialized;
  
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

  int  row_s;
  int  nnz_local, n_rows_local, n_rows_global;

  /// Index base used by specific linear solver.
  static const int index_base = 0;

  int rank, size;

  /// Private helper functions
  void
  set_lhs_structure_csr
  (int n, std::vector<int> &ir, std::vector<int> &jc);

  int
  init_lis_objects ();

  int
  assemble_lis_matrix ();
  
  int
  invoke_lis_solver ();

  void
  destroy_lis_objects ();
  
public :

 /*
 
  Typical usage :

  linear_solver *solver = new lis ();
  solver->set_lhs_structure (lhs.num_owned_rows (), ir, jc);
  solver->analyze ();
  solver->set_lhs_data (xa);
  solver->factorize ();
  solver->set_rhs (rhs);
  if (solver->solver_type () == "iterative")
  solver->set_max_iterations (1000);
  solver->set_tolerance (1e-12);
  solver->solve ();
  solver->cleanup ();


*/
  
  /// \brief Option string.
  std::string option_string;
  bool option_string_set;
  bool verbose;
  
  /// Default costructor.
  lis_distributed (LIS_INT argc = 0, char * argv[] = NULL) :
    linear_solver ("LIS", "iterative"),
    initialized (false),
    row (0),
    col (0),
    value (0),
    data (0),
    rhs (0),
    initial_guess (0),
    have_initial_guess (false),
    max_iter (1000),
    tolerance (1.0e-12),
    iterative_method ("bicg"),
    preconditioner ("none"),
    convergence_condition ("nrm2_r"),
    option_string (""),
    option_string_set (false),
    verbose (true)
  {
    lis_initialize (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
  };

  /// Set-up the matrix structure.
  void
  set_lhs_structure (int n, std::vector<int> &ir,
                     std::vector<int> &jc,
                     matrix_format_t f = csr);

  /// Perform the analysis.
  int
  analyze ();
  
  /// Set matrix entries.
  void
  set_lhs_data (std::vector<double> &xa);

  /// Set the rhs.
  void
  set_rhs (std::vector<double> &rhs_)
  { rhs = &*rhs_.begin (); }

  /// Set the initial guess.
  void
  set_initial_guess (std::vector<double> &initial_guess_)
  {
    have_initial_guess = true;
    initial_guess = &*initial_guess_.begin ();
  }

  /// Prepare the solver.
  int
  factorize ();
 

  /// Solve the system.
  int
  solve ();
 

  /// Cleanup memory.
  void
  cleanup () {
    
    destroy_lis_objects ();
    initialized = false;

  };


  /// Set maximum number of iterations (default = 1000).
  void
  set_max_iterations (int max_iter_)
  {
    max_iter = max_iter_;
    option_string_set = false;
  }

  /// Get maximum number of iterations.
  void
  get_max_iterations (int &max_iter_)
  {
    max_iter_ = max_iter;
    option_string_set = false;
  }

  /// Set tolerance of iterative method (default = 1.0e-12).
  void
  set_tolerance (double tol)
  {
    tolerance = tol;
    option_string_set = false;
  }

  /// Get tolerance of linear solver.
  void
  get_tolerance (double &tol)
  { tol = tolerance; }

  /// Get type of iterative method.
  void
  get_iterative_method (std::string &s)
  { s = iterative_method; }

  /// Get type of preconditioner.
  void
  get_preconditioner (std::string &s)
  { s = preconditioner; }

  /// get convergence condition of iterative method
  void
  get_convergence_condition (std::string &s)
  { s = convergence_condition; }



  /// Set type of iterative method (default = cg).
  void
  set_iterative_method (const std::string &s);

  /// Set type of preconditioner (default = none).
  void
  set_preconditioner (const std::string &s);

  /// Set convergence condition of iterative method
  /// (default = nrm2_r)
  void
  set_convergence_condition (const std::string &s);

  /// LIS uses 0-based indexing
  inline int
  get_index_base ()
  { return index_base; }

  
};
#endif
