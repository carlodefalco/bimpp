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
  out.resize (rows.size ());
  
  // create a set of columns to be extracted, for fast search
  // complexity: ncols * log ncols
  std::map<int, int> ordcol;
  for (jj = 0u; jj < cols.size (); ++jj)
    ordcol.insert (std::pair<int, int> (cols[jj], jj));

  for (ii = 0; ii < rows.size (); ++ii)
    if (rows[ii] < int ((*this).rows ()) && // the wanted row is actually in the matrix
        (*this)[rows[ii]].size ()) // and contains any entry
      for (auto jout = (*this)[rows[ii]].begin (); //for every column with nonzero entry in such row (constant complexity)
           jout != (*this)[rows[ii]].end ();
           ++jout)
        if (ordcol.count (jout->first)) // check if we want to extract it (log ncols)
          out[ii][ordcol.at(jout->first)] = // insert in the output matrix 
            &((*this)[rows[ii]][jout->first]);

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
  
  std::set<int> setcol;
  for (jj = 0u; jj < cols.size (); ++jj)
    setcol.insert (cols[jj]);

  for (ii = 0; ii < rows.size (); ++ii)
    if (rows[ii] < int ((*this).rows ()) &&
        (*this)[rows[ii]].size ())
      for (auto jout = (*this)[rows[ii]].begin (); 
           jout != (*this)[rows[ii]].end (); 
           ++jout)
        if (setcol.count (jout->first)) 
          out[ii][jout->first] = 
            &((*this)[rows[ii]][jout->first]);

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
