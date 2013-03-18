/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/

#ifndef HAVE_MUMPS_SOLVE_H
#define HAVE_MUMPS_SOLVE_H 1

#include <sparse.h>

void mumps_solve (sparse_matrix &lhs, std::vector<double> &rhs, const std::vector<int> &rows, const std::vector<int> &cols);

#endif
