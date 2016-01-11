/*
  Copyright (C) 2015 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file bgs_class.cpp
  \brief Block Gauss Seidel linear solver.
*/

#include "bgs_class.h"
#include <cassert>

void
bgs::init ()
{

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

}

void
bgs::set_lhs_structure
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

      blocks_size = matrix.size () / num_blocks;
      assert (blocks_size * num_blocks == matrix.size ());

      std::vector<std::vector<int> > row_indices (num_blocks);

      
      for (int ii = 0; ii < num_blocks; ++ii)
        {
          row_indices[ii].resize (blocks_size);
          for (int jj = 0; jj < blocks_size; ++jj)
            row_indices[ii][jj] = (jj + ii * blocks_size);
        }
      
      dblocks.resize (num_blocks);
      dblocks_aij.resize (num_blocks);
      ndblocks.resize (num_blocks);
  
      for (unsigned int ii = 0; ii < num_blocks; ++ii)
        {

          dblocks[ii].resize (blocks_size);
          ndblocks[ii].resize (blocks_size);
      
          std::vector<int> dcol_indices;
          std::vector<int> ndcol_indices;

          dcol_indices.resize (blocks_size);
          ndcol_indices.resize (blocks_size * (num_blocks - 1));
          auto dp = dcol_indices.begin ();
          auto ndp = ndcol_indices.begin ();
          
          for (unsigned int jj = 0; jj < matrix.size (); ++jj)
            if (jj / blocks_size == ii)
              *(dp++) = jj;
            else
              *(ndp++) = jj;

          matrix.extract_block_pointer (row_indices[ii],
                                        dcol_indices,
                                        dblocks[ii]);
      
          matrix.extract_block_pointer_keep_cols
            (row_indices[ii],
             ndcol_indices,
             ndblocks[ii]);

          dblocks[ii].aij (dblocks_aij[ii].a,
                           dblocks_aij[ii].i,
                           dblocks_aij[ii].j,
                           block_solvers[ii]->get_index_base ());
      
          block_solvers[ii]->set_lhs_structure (dblocks[ii].rows (),
                                                dblocks_aij[ii].i,
                                                dblocks_aij[ii].j);
        }
    }
  RRE = new rre (n, rre_ninit,
                 rre_nskip, rre_rank);
}

int
bgs::analyze ()
{
  for (unsigned int ii = 0; ii < num_blocks; ++ii)
    block_solvers[ii]->analyze ();
  return 0;
}

void
bgs::set_lhs_data (std::vector<double> &xa)
{
  if (rank == 0)
    {
      std::vector<double>::iterator kk = xa.begin ();
      for (sparse_matrix::row_iterator ii = matrix.begin ();
           ii != matrix.end (); ++ii)
        for (sparse_matrix::col_iterator jj = ii->begin ();
             jj != ii->end (); ++jj)
          jj->second = *(kk++);

      for (unsigned int ii = 0; ii < num_blocks; ++ii)
        {
          dblocks[ii].aij_update (dblocks_aij[ii].a,
                                  dblocks_aij[ii].i,
                                  dblocks_aij[ii].j,
                                  block_solvers[ii]->get_index_base ());
      
          block_solvers[ii]->set_lhs_data (dblocks_aij[ii].a);
        }
    }
}

void
bgs::set_rhs (std::vector<double> &rhs_)
{  
   if (rank == 0)
     {
       full_rhs = &rhs_;
       refnorm = vecnorm (full_rhs->begin (), full_rhs->end ());
       rhs.resize (num_blocks);
       for (unsigned int ii = 0; ii < num_blocks; ++ii)
         {
           rhs[ii].resize (blocks_size);
           std::copy (rhs_.begin () + ii * blocks_size,
                      rhs_.begin () + (ii + 1) * blocks_size,
                      rhs[ii].begin ());
         }
     }

}


void
bgs::set_initial_guess (std::vector<double> &guess_)
{
  have_initial_guess = true;
  if (rank == 0)
    {
      initial_guess.resize (num_blocks);
      for (unsigned int ii = 0; ii < num_blocks; ++ii)
        {
          initial_guess[ii].resize (blocks_size);
          std::copy (guess_.begin () + ii * blocks_size,
                     guess_.begin () + (ii + 1) * blocks_size,
                     initial_guess[ii].begin ());
        }
    }
}





int
bgs::factorize ()
{
  for (unsigned int ii = 0; ii < num_blocks; ++ii)
    block_solvers[ii]->factorize ();

  return 0;
}


/// Compute y = a*A*x + b*y
template<class sp>
void
sparse_dgemv  (sp& A,
              const std::vector<double>& x,
              double a, double b,
              std::vector<double> &y)
{
  typename sp::col_iterator j;
  for (unsigned int i = 0; i < A.size (); ++i)
    if (A[i].size ())
      for (j = A[i].begin (); j != A[i].end (); ++j)
        {
          y[i] *= b;
          y[i] += a * A.col_val (j) * x[A.col_idx (j)];
        }
}


int
bgs::solve ()
{
  /*
    ideally we should iteratively compute 

    x(k+1) = x(k) + P \ (b - A*x(k))

    where

    A*x(k) = a * d * x(k)

    and P is the lower block-triangular part of

    a * d

    d is a matrix composed of diagonal blocks of the form

    d = [d11    0     0    ...;
         d21    d22   0    ...;
         d31    0     d33  ...;
         ...    ...   ...  ...];

    and the computation of 

    y = P \ z

    is implemented via block-forward substitution

    y(1) = P(1,1) \ z(1) = d11 \ (a11 \ z(1))
    y(2) = P(2,2) \ (z(2) - P(2,1) * z(1)) = d22 \ (a22 \ ((a21 * d11 + a22 * d22) * z(1)))
    y(3) = P(3,3) \ (z(3) - P(3,1) * z(1) - P(3,2) * z(2)) = d33 \ (a33 \ ((a11 * d11 + a33 * d33) * z(3)))
    etc. ...

    currently we use instead

    x1(k+1) = d11 \ (a11 \ (b1 - a12 d22 x2(k) - a13 d33 x3(k)))
    x2(k+1) = d22 \ (a22 \ (b2 - (a21 d11 x1(k+1) + a22 d22 x1(k+1)) - a33 d33 x3(k))
    x3(k+1) = d33 \ (a33 \ (b3 - (a31 d11 x1(k+1) + a33 d33 x1(k+1)) - a22 d22 x2(k+1))
    etc. ...

   */
  std::vector<std::vector<double> > x(num_blocks);
  if (rank == 0)
    {

      if (! have_initial_guess)
        {
          initial_guess.resize (num_blocks);
          for (unsigned int ii = 0; ii < num_blocks; ++ii)
            {
              initial_guess[ii].resize (blocks_size);
              initial_guess[ii].assign (blocks_size, 0.0);
            }
        }
      
      for (int ii = 0; ii < num_blocks; ++ii)
        {
          x[ii].resize (blocks_size);
          std::copy (initial_guess[ii].begin (),
                     initial_guess[ii].end (),
                     full_rhs->begin () + ii * blocks_size);
        }
    }

  resnorm.resize (0);
  for (int ii = 0; ii < max_iter; ++ii)
    {
      resnorm.push_back (0);
      for (int iblock = 0; iblock < num_blocks; ++iblock)
        {
          if (rank == 0)
            {
              std::copy (rhs[iblock].begin (),
                         rhs[iblock].end (),
                         x[iblock].begin ());

              sparse_dgemv (ndblocks[iblock],
                            *full_rhs,
                            -1.0, 1.0,
                            x[iblock]);

              block_solvers[iblock]->set_rhs (x[iblock]);
            }

          block_solvers[iblock]->solve ();

          if (rank == 0)
            {
              resnorm.back () +=
                (vecdiffnorm (x[iblock].begin (),
                              x[iblock].end (),
                              full_rhs->begin () + iblock * blocks_size,
                              full_rhs->begin () + (iblock + 1) * blocks_size));
              
                std::copy (x[iblock].begin (),
                           x[iblock].end (),
                           full_rhs->begin () + iblock * blocks_size);
            }
        }

      RRE->extrapolate (*full_rhs);
      if (resnorm.back () < tolerance)
        break;
    }
  
}

void
bgs::cleanup ()
{

  for (auto ii = block_solvers.begin ();
       ii != block_solvers.end ();
       ++ii)
    (*ii)->cleanup ();

  delete RRE;
  
}
