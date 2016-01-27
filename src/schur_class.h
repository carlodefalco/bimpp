/*
  Copyright (C) 2016 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file schur_class.h
  \brief Schur complement based linear solver.
*/

#ifndef HAVE_SCHUR_CLASS_H
#define HAVE_SCHUR_CLASS_H 1

#include "bim_sparse.h"
#include "linear_solver.h"
#include <mpi.h>
#include <cstdio>
#include <fstream>

/// Schur complement linear solver.
///
/// The system is decomposed as :
///
///  [ A   a ] [ X ] = [ C ]
///  [ b   B ] [ x ] = [ c ]
///
/// and the solution is computed via
/// the following algorithm :
/// 1) g := B \ c
/// 2) G := B \ b
/// 3) A += - a * G
/// 4) C += - a * g
/// 5) X = A \ C
/// 6) x = g - G * X
///
/// Different solvers can be used
/// to invert B in steps 1) and 2)
/// than to invert A in step 5)
class schur: public linear_solver
{
  
private :

  std::vector<int> idx_A, idx_B;
  sparse_matrix matrix;
  p_sparse_matrix A, a, B, b;
  std::vector<double> C, c;
  linear_solver *solver_A, *solver_B;

  typedef struct
  {
    std::vector<double> a;
    std::vector<int> i;
    std::vector<int> j;
  } aij_struct;

  aij_struct A_aij, B_aij;
  
  int rank, size;
  static const int index_base = 0;

public :
  
  /// Init the solver instance.
  void
  init ();

  /// Default constructor.
  schur ( ) :
    linear_solver ("schur", "direct")    
  { init (); };

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

  void
  set_initial_guess (std::vector<double> &guess_);

  /// Perform the factorization.
  int
  factorize ();

  /// Perform the back-substitution.
  int
  solve ();

  /// Cleanup memory.
  void
  cleanup ();

  /// uses 0-based indexing
  inline int
  get_index_base ()
  { return index_base; }

  
  /// print block-decomposed matrix
  void
  print_blocks (std::string basename = "block");
  
};
  
#endif
