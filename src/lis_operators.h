/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/


#ifndef HAVE_LIS_OPERATOR_H
#define HAVE_LIS_OPERATOR_H 1

#include <lis.h>
#include <bim_sparse.h>

struct linear_solver_option
{
  char* tolerance;
  char* maxit;
  char* other_opt;
};

void
lis_matrix_parallelization(sparse_matrix& sp, 
			   sparse_matrix& sp_loc);

void
lis_vector_parallelization(std::vector<double>& v,
			   std::vector<double>& v_loc);

void
lis_vector_unification(LIS_VECTOR& v_loc,
		       std::vector<double>& v,
		       int vsize,
		       int istart,	
		       int iend);

void
lis_solve_system(sparse_matrix& lhs, 
		 std::vector<double>& rhs, 
		 std::vector<double>& sol,
		 LIS_INT& iter, 
		 double& time, 
		 int nnodes,
		 linear_solver_option& option);
#endif
