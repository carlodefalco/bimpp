/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file lis_distributed_class.cpp
  \brief interface for linear solver built with for lis library.
*/

#include "lis_distributed_class.h"
#include <stdlib.h>
#include <sstream>
#include <string>
#include <cstring>

void
lis_distributed::set_lhs_structure_aij
(int n,
 std::vector<int> &ir,
 std::vector<int> &jc)
{
  n_row = n;
  nnz = ir.size ();

  //row = (LIS_INT *) malloc (nnz * sizeof(LIS_INT));
  //col = (LIS_INT *) malloc (nnz * sizeof(LIS_INT));
  //value = (LIS_SCALAR *) malloc (nnz * sizeof(LIS_SCALAR));
  
  row = new LIS_INT[nnz];
  col = new LIS_INT[nnz];
  value = new LIS_SCALAR[nnz];

 // lis_matrix_create (MPI_COMM_WORLD, &A);
  //lis_matrix_set_size (A, n_row, 0);  // sull'ultimo non sono sicura, forse è n_row
  for (int i = 0; i < nnz; i++)
    {
       row[i] = ir[i];
       col[i] = jc[i] - index_base;
       value[i] = 0.0;
    }
  //lis_matrix_set_coo (nnz, row, col, value, A);
}

void
lis_distributed::set_lhs_data (std::vector<double> &xa)
{
  data = &*xa.begin ();
}

void
lis_distributed::cleanup_master () { }


