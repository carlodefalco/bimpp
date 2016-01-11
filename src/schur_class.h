/*
  Copyright (C) 2015 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file schur_class.h
  \brief wrapper for mumps data.
*/

#ifndef HAVE_SCHUR_CLASS_H
#define HAVE_SCHUR_CLASS_H 1

#include "bim_sparse.h"
#include "linear_solver.h"

/// Class that solves a linear system using an intermediate Schur-complement elimination step.
class schur: public linear_solver
{
private :

  linear_solver* master;
  linear_solver* slave;
  
public :

  /// Init the (serial) mumps solver instance.
  void
  init ();

  /// Default constructor.
  schur (bool verbose_ = false, int icntl23_ = 0, int working_host_ = 1) :
    linear_solver ("MUMPS", "direct"),
    verbose (verbose_),
    icntl23 (icntl23_),
    working_host (working_host_)
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

  /// Perform the factorization.
  int
  factorize ();

  /// Perform the back-substitution.
  int
  solve ();

  /// Cleanup memory.
  void
  cleanup ();

  /// MUMPS uses 1-based indexing
  inline int
  get_index_base ()
  { return index_base; }

};

#endif
