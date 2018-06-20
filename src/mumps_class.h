/*
 *  Copyright (C) 2011, 2015 Carlo de Falco
 *  This software is distributed under the terms
 *  the terms of the GNU/GPL licence v3
 */
/*! \file mumps_class.h
 *  \brief wrapper for mumps data.
 */

#ifndef HAVE_MUMPS_CLASS_H
#define HAVE_MUMPS_CLASS_H 1

#define F77_COMM_WORLD -987654
#define JOB_INIT -1
#define JOB_ANALYZE 1
#define JOB_FACTORIZE 2
#define JOB_SOLVE 3
#define JOB_END  -2

#include "bim_sparse.h"
#include "bim_distributed_vector.h"
#include <dmumps_c.h>
#include "linear_solver.h"

#include <memory>

//using namespace bim;

/// Wrapper class around the MUMPS linear solver.
class mumps: public linear_solver
{
private :
  
  bool verbose;
  int  icntl23;
  int working_host;
  static const int index_base = 1;
  
  std::shared_ptr<distributed_vector> glob_rhs;
  
public :
  
  DMUMPS_STRUC_C id;
  
  /// Init the mumps solver instance.
  void
  init ();
  
  /// Default constructor.
  mumps (bool verbose_ = false, int icntl23_ = 0, int working_host_ = 1) :
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
  
  /// Set-up the structure of a distributed matrix.
  void
  set_distributed_lhs_structure
  (int n,
   std::vector<int> &ir,
   std::vector<int> &jc,
   matrix_format_t f = aij);
  
  /// Set the matrix format to distributed.
  void
  set_lhs_distributed ();
  
  /// Perform the analysis.
  int
  analyze ();
  
  /// Set matrix entries.
  void
  set_lhs_data (std::vector<double> &xa);
  
  /// Set matrix entries of a distributed matrix.
  void
  set_distributed_lhs_data (std::vector<double> &xa);
  
  /// Set the rhs.
  void
  set_rhs (std::vector<double> &rhs);
  
  /// Set the rhs from a distributed_vector.
  void
  set_rhs_distributed (distributed_vector &rhs);
  
  /// Perform the factorization.
  int
  factorize ();
  
  /// Perform the back-substitution.
  int
  solve ();

  /// Returns the distributed_vector solution.
  distributed_vector &
  get_distributed_solution ()
  {
    assert (glob_rhs != nullptr);
    
    // Communicate non-local values to their owners first.
    glob_rhs->assemble (replace_op);
    
    return *glob_rhs;
  }
  
  /// Cleanup memory.
  void
  cleanup ();
  
  /// MUMPS uses 1-based indexing
  inline int
  get_index_base ()
  { return index_base; }
  
};

#endif
