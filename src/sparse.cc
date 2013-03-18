/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/

#include <sparse.h>

template<>
double double_sparse_matrix::col_val (double_sparse_matrix::col_iterator j)   {return (*j).second;}

template<>
double double_p_sparse_matrix::col_val (double_p_sparse_matrix::col_iterator j) {return *((*j).second);}

void sparse_matrix::extract_block_pointer (const std::vector<int> &rows, const std::vector<int> &cols, p_sparse_matrix &out)
{
  int  ii, jj;
  out.resize (rows.size ());
  
  for (ii = 0; ii < rows.size (); ++ii)    
    if (rows[ii] < (*this).rows () && (*this)[rows[ii]].size ())      
      for (jj = 0; jj < cols.size (); ++jj)
        if ((*this)[rows[ii]].count (cols[jj]))
          out[ii][jj] = &((*this)[rows[ii]][cols[jj]]);

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
