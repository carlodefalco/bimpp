/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file mumps_solve.cpp
  \brief wrapper function for a complete linear system solution.
*/

#include <mumps_solve.h>
#include <mumps_class.h>

void
mumps_solve (sparse_matrix &lhs,
             std::vector<double> &rhs,
             const std::vector<int> &rows,
             const std::vector<int> &cols)
{

  std::vector<double> lrhs;
  lrhs.resize (cols.size ());

  for (size_t i = 0; i < cols.size (); ++i)
    lrhs[i] = rhs[cols[i]];

  std::vector<int> ir, jc;
  std::vector<double> xa;

  p_sparse_matrix llhs;
  lhs.extract_block_pointer (rows, cols, llhs);
  llhs.aij (xa, ir, jc, 1);

  mumps mumps_solver;

  mumps_solver.set_lhs_structure (llhs.rows (), ir, jc);
  mumps_solver.analyze ();

  mumps_solver.set_lhs_data (xa);
  mumps_solver.factorize ();

  mumps_solver.set_rhs (lrhs);
  mumps_solver.solve ();

  mumps_solver.cleanup ();

  for (size_t i = 0; i < lrhs.size (); ++i)
    rhs[cols[i]] = lrhs[i];

}
