/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include "bim_sparse.h"

//using namespace bim;

/// Template specialization for double matrix
// template<>
// double
// double_sparse_matrix::col_val (double_sparse_matrix::col_iterator j)   {return (*j).second;}

/// Template specialization for double* matrix
// template<>
// double
// double_p_sparse_matrix::col_val (double_p_sparse_matrix::col_iterator j) {return *((*j).second);}


void
sparse_matrix::extract_block_pointer (const std::vector<int> &rows,
                                      const std::vector<int> &cols,
                                      p_sparse_matrix &out)
{
  size_t  ii, jj;
  int nrows = rows.size ();
  out.resize (nrows);

  for (ii = 0; ii < (nrows); ++ii)
    if (rows[ii] < int ((*this).rows ()) &&
        (*this)[rows[ii]].size ())
      for (jj = 0; jj < cols.size (); ++jj)
        if ((*this)[rows[ii]].count (cols[jj]))
          out[ii][jj] = & ((*this)[rows[ii]][cols[jj]]);

  out.set_properties ();
}

void
sparse_matrix::extract_block_pointer_keep_cols
(const std::vector<int> &rows,
 const std::vector<int> &cols,
 p_sparse_matrix &out)
{
  size_t  ii, jj;
  out.resize (rows.size ());

  for (ii = 0; ii < rows.size (); ++ii)
    if (rows[ii] < int ((*this).rows ()) &&
        (*this)[rows[ii]].size ())
      for (jj = 0; jj < cols.size (); ++jj)
        if ((*this)[rows[ii]].count (cols[jj]))
          out[ii][cols[jj]] = & ((*this)[rows[ii]][cols[jj]]);

  out.set_properties ();
}

void sparse_matrix::reset ()
{
  double_sparse_matrix::row_iterator ii;
  double_sparse_matrix::col_iterator jj;
  for (ii = this->begin (); ii != this->end (); ++ii)
    for (jj = (*ii).begin (); jj != (*ii).end (); ++jj)
      (*jj).second = 0.0;
}

std::vector<double>
operator * (sparse_matrix& M, const std::vector<double>& x)
{
  std::vector<double> y (M.rows (), 0.0);
  sparse_matrix::col_iterator j;
  for (unsigned int i = 0; i < M.size (); ++i)
    if (M[i].size ())
      for (j = M[i].begin (); j != M[i].end (); ++j)
        y[i] += M.col_val (j) * x[M.col_idx (j)];

  return y;
}
