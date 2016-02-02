/*
  Copyright (C) 2016 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file schur_class.cpp
  \brief Schur complement based linear solver.
*/

#include "schur_class.h"

void
schur::init ()
{
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
}


void
schur::set_lhs_structure
(int n,
 std::vector<int> &ir,
 std::vector<int> &jc,
 matrix_format_t f)
{

  if (rank == 0)
    {
      matrix.resize (n);
      if (f == csr)
        for (unsigned int i = 0; i < n; ++i)
          for (unsigned int j = ir[i]; j < ir[i+1]; ++j)
            matrix[i][jc[j]] = 0.0;
      else
        {
          assert (ir.size () == jc.size ());
          for (unsigned int i = 0; i < ir.size (); ++i)
            matrix[ir[i]][jc[i]] = 0.0;
        }
      
      matrix.extract_block_pointer (idx_A, idx_A, A);
      matrix.extract_block_pointer (idx_A, idx_B, a);
      matrix.extract_block_pointer (idx_B, idx_A, b);
      matrix.extract_block_pointer (idx_B, idx_B, B);
      
      A.aij (A_aij.a, A_aij.i, A_aij.j,
             solver_A->get_index_base ());
      
      B.aij (B_aij.a, B_aij.i, B_aij.j,
             solver_B->get_index_base ());
      
    }
}

int
schur::analyze ()
{
  solver_A->analyze ();
  solver_B->analyze ();
  return 0;
}

void
schur::set_lhs_data (std::vector<double> &xa)
{
  if (rank == 0)
    {
      // update aij format
      auto kk = xa.begin ();
      for (auto rr = matrix.begin (); rr != matrix.end (); ++rr)
        for (auto cc = rr->begin (); cc != rr->end (); ++cc)
          cc->second = *(kk++);
      
      // update diagonal blocks and solvers
      A.aij_update (A_aij.a, A_aij.i, A_aij.j,
                    solver_A->get_index_base ());
      solver_A->set_lhs_data (A_aij.a);

      B.aij_update (B_aij.a, B_aij.i, B_aij.j,
                    solver_B->get_index_base ());
      solver_B->set_lhs_data (B_aij.a);
    }
  return;
}

void
schur::set_rhs (std::vector<double> &rhs)
{
  if (rank == 0)
    {
      C.resize (idx_A.size ());
      auto idx_C = C.begin ();
      for (int ii = 0; ii < idx_A.size (); ++ii)
        *(idx_C++) = rhs[idx_A[ii]];
        
      c.resize (idx_B.size ());
      auto idx_c = c.begin ();
      for (int ii = 0; ii < idx_B.size (); ++ii)
        *(idx_c++) = rhs[idx_B[ii]];     
    }

  return;
}

int
schur::factorize ()
{
  int res_A = solver_A->factorize ();
  int res_B = solver_B->factorize ();
  return (res_A == 0 && res_B == 0) ? 0 : -1;
}

int
schur::solve () { return 0; }

void
schur::cleanup () { return; }
