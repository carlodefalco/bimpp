/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file lis_class.h
  \wrapper for lis data.
*/

#ifndef HAVE_LIS_CLASS
#define HAVE_LIS_CLASS 1

#include <lis.h>
#include <bim_sparse.h>
#include <string>
#include <linear_solver.h>

class lis: public linear_solver
{
private :
  std::vector<int> row_ptr;
  std::vector<int> jcol;
  std::vector<double> data;
  std::vector<double> rhs;
  std::vector<double>::iterator rhs_it;

  int max_iter;
  double tolerance;

  std::string lin_solver;
  std::string preconditioner;
  std::string other_options;

  int i_s, row_s;
  int n, nnz, n_row;

public :
  /// Default costructor.
  lis (LIS_INT argc = 0, char * argv[] = NULL) :
    linear_solver ("LIS", "iterative")
  {
    lis_initialize (&argc, &argv);
    max_iter = 1000;
    tolerance = 1.0e-12;
  };

  /// Set-up the matrix structure.
  void
  set_lhs_structure (int n,
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

  /// Perform the factorization.
  int
  factorize ()
  {
    return 1;
  };

  /// Solve the system.
  int
  solve ();

  /// Cleanup memory.
  void
  cleanup ();

  /// Set maximum number of iterations.
  void
  set_max_iterations (int max_iter_);

  /// Get maximum number of iterations.
  void
  get_max_iterations (int &max_iter_);

  /// Set tolerance of linear solver.
  void
  set_tolerance (double tol);

  /// Get tolerance of linear solver.
  void
  get_tolerance (double &tol);

  /// Set type of linear solver.
  void
  set_linear_solver (const std::string &s);

  /// Get type of linear solver.
  void
  get_linear_solver (std::string &s);

  /// Set type of preconditioner.
  void
  set_preconditioner (const std::string &s);

  /// Get type of preconditioner.
  void
  get_preconditioner (std::string &s);

  /// Set other optins of linear solver
  void
  set_other_options (const std::string &s);

  /// get other options of linear solver
  void
  get_other_options (std::string &s);

};
#endif
