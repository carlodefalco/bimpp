/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/
/*! \file mumps_class.h
  \brief wrapper for mumps data.
*/

#ifndef HAVE_MUMPS_CLASS_H
#define HAVE_MUMPS_CLASS_H 1

#define F77_COMM_WORLD -987654
#define JOB_INIT -1
#define JOB_ANALYZE 1
#define JOB_FACTORIZE 2
#define JOB_SOLVE 3
#define JOB_END  -2

#include <bim_sparse.h>
#include <dmumps_c.h>

//using namespace bim;

/// Wrapper class around the MUMPS linear solver.
class mumps
{
private :
  bool verbose;
  int  icntl23;

public :
  DMUMPS_STRUC_C id;
  
  /// Init the (serial) mumps solver instance.
  void 
  init ();

  /// Default constructor.
  mumps (bool verbose_ = false, int icntl23_ = 0) : verbose (verbose_), icntl23 (icntl23_) 
  {init ();};

  /// Set-up the matrix structure.
  void 
  set_lhs_structure (int n, 
                     std::vector<int> &ir, 
                     std::vector<int> &jc);

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
};



#endif
