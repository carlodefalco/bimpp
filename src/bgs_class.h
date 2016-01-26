/*
  Copyright (C) 2015 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file bgs_class.h
  \brief Block Gauss Seidel linear solver.
*/

#ifndef HAVE_BGS_CLASS_H
#define HAVE_BGS_CLASS_H 1

#include "bim_sparse.h"
#include "linear_solver.h"
#include "rre.h"
#include <mpi.h>
#include <cstdio>
#include <fstream>
#include <iostream>

/// Block Gauss Seidel linear solver.
class bgs: public linear_solver
{
private :

  /// Number of blocks.
  const unsigned int num_blocks;
  /// Size of blocks.
  unsigned int blocks_size;
  
  /// Solvers for each block.
  std::vector<linear_solver*> block_solvers;
  /// Max iterations in solve [default = 1000].
  int max_iter;
  /// Applied relative tolerance [default = 1e-12].
  double tolerance;

  /// Full matrix.
  sparse_matrix matrix;
  /// Diagonal blocks.
  std::vector<p_sparse_matrix> dblocks;
  /// Nondiagonal blocks.
  std::vector<p_sparse_matrix> ndblocks;
  /// Source term blocks.
  std::vector<std::vector<double> > rhs;
  /// Full source term.
  std::vector<double> *full_rhs;
  /// Initial guess blocks. 
  std::vector<std::vector<double> > initial_guess;
  
  /// Struct for aij format.
  typedef struct
  {
    std::vector<double> a;
    std::vector<int> i;
    std::vector<int> j;
  } aij_struct;

  /// Diagonal blocks converted to aij.
  std::vector<aij_struct> dblocks_aij;

  /// Flag for getting initial guess.
  bool have_initial_guess;
  
  /// MPI rank.
  int rank; 
  /// MPI size.
  int size;

  /// History of residual norms.
  std::vector<double> resnorm;
  /// Reference value.
  double refnorm;
  /// Parameters for rre.
  int rre_ninit, rre_nskip, rre_rank;
  /// Instance of rre extrapolator.
  rre *RRE;
  static const int index_base = 0;

  /// Values for right preconditioner
  std::vector<double> rprec;

  /// Storage for row indices
  std::vector<int> ir;
  /// Storage for column indices
  std::vector<int> jc;
  /// Storage for structure type (csr or aij)
  matrix_format_t matf;

public :
  
  /// Init the solver instance.
  void
  init ();

  /// Default constructor.
  /**
   The instance of bgs accepts a vector of
   linear solvers, which defines the number of blocks
   in which the system will be split. Optional parameters
   are used in setting up the RR extrapolation.
   */
  bgs (const std::vector<linear_solver*> &block_solvers_,
       int ninit_ = 2, int nskip_ = 0, int rank_ = 3) :
    linear_solver ("BGS", "iterative"),
    num_blocks (block_solvers_.size ()),
    block_solvers (block_solvers_),
    max_iter (1000),
    tolerance (1.0e-12),
    have_initial_guess (false),
    rre_ninit (ninit_),
    rre_nskip (nskip_),
    rre_rank (rank_),
    RRE (0)
  { init (); };

  /// Set-up the matrix structure.
  /**
   The incoming matrix structure is divided in uniform 
   blocks. The number of blocks is given by the number 
   of solvers loaded by the constructor.
   */
  void
  set_lhs_structure
  (int n,
   std::vector<int> &ir,
   std::vector<int> &jc,
   matrix_format_t f = aij);

  /// Perform the analysis.
  int
  analyze ();

  /// Set right preconditioner data.
  /**
   Set the right preconditioner data, which 
   is supposed to be:
   - made of diagonal blocks
   - have nonzero blocks on the diagonal and first column only
   and therefore is represented as a vector, with the diagonals
   stored successively (d11, d22, d33 ... d21, d31 ...)
   */
  void
  set_preconditioner_data (std::vector<double> &xd);

  /// Set matrix entries.
  void
  set_lhs_data (std::vector<double> &xa);

  /// Set the rhs.
  void
  set_rhs (std::vector<double> &rhs);

  /// Set the initial guess.
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

  /// Uses 0-based indexing
  inline int
  get_index_base ()
  { return index_base; }

  
  /// Print block-decomposed matrix
  void
  print_blocks (std::string basename = "block")
  {
    for (unsigned int ii = 0; ii < num_blocks; ++ii)
      {
        char tmp[255];

        sprintf (tmp, "%s_%d_diagonal.m", basename.c_str (), ii);
        std::ofstream fout (tmp);
        fout << dblocks[ii];
        fout.close ();

        sprintf (tmp, "%s_%d_non_diagonal.m", basename.c_str (), ii);
        fout.open (tmp);
        fout << ndblocks[ii];
        fout.close ();

      }
  }
  
  /// Set maximum number of iterations (default = 1000).
  void
  set_max_iterations (int max_iter_)
  { max_iter = max_iter_; }
  
  /// Get maximum number of iterations.
  void
  get_max_iterations (int &max_iter_)
  { max_iter_ = max_iter; }

  /// Set tolerance of iterative method (default = 1.0e-12).
  void
  set_tolerance (double tol)
  { tolerance = tol; }
  
  /// Get tolerance of linear solver.
  void
  get_tolerance (double &tol)
  { tol = tolerance; }

  /// Utility for computing vector norms.
  double
  vecnorm (std::vector<double>::iterator first,
           std::vector<double>::iterator last)
  {
    double n = 0;
    for (auto i = first; i != last; ++i)
      n += (*i) * (*i);
    return (n);
  }

  /// Utility for computing vector distances.
  double
  vecdiffnorm (std::vector<double>::iterator xfirst,
               std::vector<double>::iterator xlast,
               std::vector<double>::iterator yfirst,
               std::vector<double>::iterator ylast)
  {
    double n = 0, tt = 0;
    for (auto i = xfirst, j = yfirst;
         i != xlast && j != ylast;
         ++i, ++j)
      {
        tt = (*i) - (*j);
        n += tt * tt;
      }
    return (n);
  }


  /// Print residual norm history.
  void
  print_resnorm (std::string basename = "resnorm")
  {
    char tmp[255];
    
    sprintf (tmp, "%s.m", basename.c_str ());
    std::ofstream fout (tmp);

    int ii = 1;
    for (auto i = resnorm.begin (); i != resnorm.end (); ++i)
      fout << "r(" << ii++ << ") = " << (*i) << ";" << std::endl;

    fout.close ();
    
  }
  
};

#endif
