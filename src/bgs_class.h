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

/// Block Gauss Seidel linear solver.
class bgs: public linear_solver
{
private :

  const unsigned int num_blocks;
  unsigned int blocks_size;
  
  std::vector<linear_solver*> block_solvers;
  // [default = 1000].
  int max_iter;
  // [default = 1e-12].
  double tolerance;

  sparse_matrix matrix;
  std::vector<p_sparse_matrix> ndblocks, dblocks;
  std::vector<std::vector<double> > rhs;
  std::vector<double> *full_rhs;
  std::vector<std::vector<double> > initial_guess;
  
  typedef struct
  {
    std::vector<double> a;
    std::vector<int> i;
    std::vector<int> j;
  } aij_struct;

  std::vector<aij_struct> dblocks_aij;

  bool have_initial_guess;
  int rank, size;

  std::vector<double> resnorm;
  double refnorm;
  int rre_ninit, rre_nskip, rre_rank;
  rre *RRE;
  static const int index_base = 0;

public :
  
  /// Init the solver instance.
  void
  init ();

  /// Default constructor.
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
    rre_rank (rank_)
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
  print_blocks (std::string basename = "block")
  {
    for (int ii = 0; ii < num_blocks; ++ii)
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

  double
  vecnorm (std::vector<double>::iterator first,
           std::vector<double>::iterator last)
  {
    double n = 0;
    for (auto i = first; i != last; ++ i)
      n += (*i) * (*i);
    return (n);
  }

  double
  vecdiffnorm (std::vector<double>::iterator xfirst,
               std::vector<double>::iterator xlast,
               std::vector<double>::iterator yfirst,
               std::vector<double>::iterator ylast)
  {
    double n = 0, tt = 0;
    for (auto i = xfirst, j = yfirst;
         i != xlast && j != ylast;
         ++ i, ++ j)
      {
        tt = (*i) - (*j);
        n += tt * tt;
      }
    return (n);
  }


  /// print block-decomposed matrix
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
