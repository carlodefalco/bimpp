/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/

#ifndef HAVE_MUMPS_CLASS_H
#define HAVE_MUMPS_CLASS_H 1

#define F77_COMM_WORLD -987654
#define JOB_INIT -1
#define JOB_ANALYZE 1
#define JOB_FACTORIZE 2
#define JOB_SOLVE 3
#define JOB_END  -2

#include <sparse.h>
#include <dmumps_c.h>

/// Wrapper class around the MUMPS linear solver.
class mumps
{
 public :
  DMUMPS_STRUC_C id;
  
  /// Init the (serial) mumps solver instance.
  void init ();

  /// Default constructor.
  mumps () {init ();};

  /// Set-up the matrix structure.
  void set_lhs_structure (int n, std::vector<int> &ir, std::vector<int> &jc);

  /// Perform the analysis.
  void analyze ();

  /// Set matrix entries.
  void set_lhs_data (std::vector<double> &xa);

  /// Set the rhs.
  void set_rhs (std::vector<double> &rhs);

  /// Perform the factorization.
  void factorize ();

  /// Perform the back-substitution.
  void solve ();

  /// Cleanup memory.
  void cleanup ();
};



#endif
