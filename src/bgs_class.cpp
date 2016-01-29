/*
  Copyright (C) 2015 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file bgs_class.cpp
  \brief Block Gauss Seidel linear solver.
*/

#include "bgs_class.h"
#include "bim_timing.h"
#include <cassert>
#include <iostream>

void
bgs::init ()
{

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

}

void
bgs::set_lhs_structure
(int n,
 std::vector<int> &_ir,
 std::vector<int> &_jc,
 matrix_format_t f)
{
  if (rank == 0)
    {
      // store original structure
      ir.resize (_ir.size ()); std::copy (_ir.begin (), _ir.end (), ir.begin ());
      jc.resize (_jc.size ()); std::copy (_jc.begin (), _jc.end (), jc.begin ());
      matf = f;

      // compute block size
      matrix.resize (n);

      blocks_size = matrix.size() / num_blocks;
      assert (blocks_size * num_blocks == matrix.size ());

      // create the sparse matrix
      if (f == csr)
        for (unsigned int i = 0; i < (unsigned int)n; ++i)
          for (unsigned int j = ir[i]; 
               j < (unsigned int)(ir[i+1]); ++j)
          {
            matrix[i][jc[j]] = 0.0;
            // add structure due to preconditioning
            if (jc[j] - blocks_size >= 0)
              matrix[i][jc[j] % blocks_size] = 0.0;
          }
      else
        {
          assert (ir.size () == jc.size ());
          for (unsigned int i = 0; i < ir.size (); ++i)
          {
            matrix[ir[i]][jc[i]] = 0.0;
            // add structure due to preconditioning
            if (jc[i] - blocks_size >= 0)
              matrix[ir[i]][jc[i] % blocks_size] = 0.0;
          }
        }
      
      // create map from [block,localrow] to [globalrow]
      std::vector<std::vector<int> > row_indices (num_blocks);
      for (unsigned int ii = 0; ii < num_blocks; ++ii)
        {
          row_indices[ii].resize (blocks_size);
          for (unsigned int jj = 0; jj < blocks_size; ++jj)
            row_indices[ii][jj] = (jj + ii * blocks_size);
        }
      
      // create blocks
      dblocks.resize (num_blocks);
      dblocks_aij.resize (num_blocks);
      ndblocks.resize (num_blocks);

      // create maps from local matrices to full matrix
      for (unsigned int ii = 0; ii < num_blocks; ++ii)
        {

      tic ();
          // set num rows
          dblocks[ii].resize (blocks_size);
          ndblocks[ii].resize (blocks_size);
      toc ("set num rows");

      tic ();
          // get indices
          std::vector<int> dcol_indices;
          std::vector<int> ndcol_indices;
      toc ("get indices");

      tic ();
          dcol_indices.resize (blocks_size);
          ndcol_indices.resize (blocks_size * (num_blocks - 1));
          auto dp = dcol_indices.begin ();
          auto ndp = ndcol_indices.begin ();
      toc ("init vars");

      tic ();
          // if on diagonal columns for this block, 
          // add column in diagonal map
          for (unsigned int jj = 0; jj < matrix.size (); ++jj)
            if (jj / blocks_size == ii)
              *(dp++) = jj;
            else
              *(ndp++) = jj;
          // else add column in nondiagonal map
      toc ("split columns");

      tic ();
          //extract diagonal block as a square matrix
          matrix.extract_block_pointer (row_indices[ii],
                                        dcol_indices,
                                        dblocks[ii]);
      toc ("extract diagonal");

      tic ();
          // extract nondiagonal block as blocks_size * n matrix
          matrix.extract_block_pointer_keep_cols
            (row_indices[ii],
             ndcol_indices,
             ndblocks[ii]);
      toc ("extract nondiagonal");

          // convert diagonal blocks to aij
      tic ();
          dblocks[ii].aij (dblocks_aij[ii].a,
                           dblocks_aij[ii].i,
                           dblocks_aij[ii].j,
                           block_solvers[ii]->get_index_base ());
      toc ("convert diag blocks");

          // set diagonal block internal structure 
      tic ();
          block_solvers[ii]->set_lhs_structure (dblocks[ii].rows (),
                                                dblocks_aij[ii].i,
                                                dblocks_aij[ii].j);
      toc ("pass blocks to resp. solver");
        }

      // init rre 
      RRE = new rre (n, rre_ninit,
                     rre_nskip, rre_rank);

      rprec.assign ((2 * num_blocks - 1) * blocks_size, 1.0);
      auto ite = rprec.begin () + (num_blocks * blocks_size);
      std::fill (ite, rprec.end (), 0.0);
       
    }
}

int
bgs::analyze ()
{
  for (unsigned int ii = 0; ii < num_blocks; ++ii)
    block_solvers[ii]->analyze ();
  return 0;
}

void
bgs::set_preconditioner_data (std::vector<double> &xd)
{
  if (rank == 0)
    {
      assert (xd.size () == ((2 * num_blocks - 1) * blocks_size));
      std::copy (xd.begin (), xd.end (), rprec.begin ());
    }
}

void
bgs::set_lhs_data (std::vector<double> &xa)
{
  if (rank == 0)
    {
      assert (xa.size () == jc.size ());

      // cleanup data needed if original and
      // preconditioned matrix have different
      // structure
      for (auto rr = matrix.begin (); rr != matrix.end (); ++rr)
        for (auto cc = rr->begin (); cc != rr->end (); ++cc)
            cc->second = 0.0;
      
      // update the sparse matrix
      if (matf == csr)
        for (unsigned int i = 0; i < matrix.size (); ++i)
          for (unsigned int j = ir[i]; 
               j < (unsigned int)(ir[i+1]); ++j)
          {
            // diagonal scaling
            matrix[i].at (jc[j]) = xa[j] * rprec[jc[j]];
            // nondiagonal part of precond.
            if (jc[j] >= (int) blocks_size)
              matrix[i].at (jc[j] % blocks_size) += 
                xa[j] * rprec[jc[j] + 
                  blocks_size * (num_blocks - 1)];
          }
      else // if matf == aij
        {
          for (unsigned int i = 0; i < ir.size (); ++i)
          {
            // diagonal scaling
            matrix[ir[i]].at (jc[i]) = xa[i] * rprec[jc[i]];
            // nondiagonal part of precond.
            if (jc[i] >= (int) blocks_size)
              matrix[ir[i]].at (jc[i] % blocks_size) += 
                xa[i] * rprec[jc[i] + 
                  blocks_size * (num_blocks - 1)];
          }
        }

      // update diagonal blocks and solvers
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

  // /* The inverse of the preconditioner:
  //  * 
  //  * d11 0   0   0   ...
  //  * d21 d22 0   0   ...
  //  * d31 0   d33 0   ...
  //  * d41 0   0   d44 ...
  //  * ...
  //  *
  //  * is
  //  *
  //  * 1/d11         0     0     0     ...
  //  * -d21/(d22d11) 1/d22 0     0     ...
  //  * -d31/(d33d11) 0     1/d33 0     ...
  //  * -d41/(d44d11) 0     0     1/d44 ... 
  //  */
  // auto xx = full_rhs->begin ();
  // // restore the first block
  // for (unsigned int kk = 0; kk < blocks_size; kk++, xx++)
  //   *xx = *xx / rprec[kk];
  // // restore further blocks
  // auto precd = rprec.begin () + blocks_size;
  // auto precnd = rprec.begin () + num_blocks * blocks_size;
  // for (unsigned int ii = 1; ii < num_blocks; ii++)
  //   {
  //     auto zz = full_rhs->begin ();
  //     for (unsigned int jj = 0; jj < blocks_size; 
  //       jj++, xx++, zz++, precd++, precnd++)
  //       *xx = (*xx - *precnd * *zz) / (*precd);
  //   }

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

    d is a block matrix of the form

    d = [d11    0     0    ...;
         d21    d22   0    ...;
         d31    0     d33  ...;
         ...    ...   ...  ...];

    where blocks are all diagonal matrices; 
    the computation of 

    y = P \ z

    is implemented via block-forward substitution

    y(1) = P(1,1) \ z(1) 
    y(2) = P(2,2) \ (z(2) - P(2,1) * z(1))
    y(3) = P(3,3) \ (z(3) - P(3,1) * z(1) - P(3,2) * z(2)) 
    etc. ...

    currently we use instead

    x1(k+1) = (A11 \ (b1 - A12 x2(k)   - A13 x3(k)))
    x2(k+1) = (A22 \ (b2 - A21 x1(k+1) - A33 x3(k)))
    x3(k+1) = (A33 \ (b3 - A31 x1(k+1) - A22 x2(k+1)))
    etc. ...

   */

  std::vector<std::vector<double> > x(num_blocks);
  int retval = 0;
  int stop_bgs_loop = 0;
    
  if (rank == 0)
    {

      // create initial guess if needed
      if (! have_initial_guess)
        {
          initial_guess.resize (num_blocks);
          for (unsigned int ii = 0; ii < num_blocks; ++ii)
            {
              initial_guess[ii].resize (blocks_size);
              initial_guess[ii].assign (blocks_size, 0.0);
            }
        }

      // set initial guess
      for (unsigned int ii = 0; ii < num_blocks; ++ii)
        {
          x[ii].resize (blocks_size);
          std::copy (initial_guess[ii].begin (),
                     initial_guess[ii].end (),
                     full_rhs->begin () + ii * blocks_size);
        }
    }

  resnorm.resize (0);
  for (int iter = 0; iter < max_iter; ++iter)
    {
      resnorm.push_back (0);
      for (unsigned int ii = 0; ii < num_blocks; ++ii)
        {
          if (rank == 0)
            {
              // copy rhs onto x
              std::copy (rhs[ii].begin (),
                         rhs[ii].end (),
                         x[ii].begin ());

              // x[ii] <- -nd[ii]*fullrhs + x[ii]
              sparse_dgemv (ndblocks[ii],
                            *full_rhs,
                            -1.0, 1.0,
                            x[ii]);

              block_solvers[ii]->set_rhs (x[ii]);
            }

          retval = block_solvers[ii]->solve ();

          if (rank == 0)
            {
              resnorm.back () +=
                (vecdiffnorm (x[ii].begin (),
                              x[ii].end (),
                              full_rhs->begin () + ii * blocks_size,
                              full_rhs->begin () + 
                                (ii + 1) * blocks_size));
              
                std::copy (x[ii].begin (),
                           x[ii].end (),
                           full_rhs->begin () + ii * blocks_size);
               
            }
        }


      if (rank == 0)
        {
          if (resnorm.back () < tolerance)
            stop_bgs_loop = 1;
          else
            RRE->extrapolate (*full_rhs);                  
        }

      MPI_Bcast (&stop_bgs_loop, 1, MPI_INT,
                 0, MPI_COMM_WORLD);

     if (stop_bgs_loop == 1)
       break;
    }
  
  /*
   * GIVEN THE SHAPE OF PRECONDITIONER,
   * recovering the original variables consists in
   * scaling with the preconditioner:
   *
   * (A) the first block is given by
   *           x1 = d11 * y1
   * (B) the remaining blocks, by
   *           xk = dk1 * y1 + dkk yk
   *  
   * in order to work directly on the array, 
   * step (B) is performed first, step (A) later.
   *   
   */
   auto xx = full_rhs->begin () + blocks_size;
   // restore latter blocks
   auto precd = rprec.begin () + blocks_size;
   auto precnd = rprec.begin () + num_blocks * blocks_size;
   for (unsigned int ii = 1; ii < num_blocks; ii++)
     {
       auto zz = full_rhs->begin ();
       for (unsigned int jj = 0; jj < blocks_size; 
         jj++, xx++, zz++, precd++, precnd++)
         *xx = *precd * *xx + *precnd * *zz;
     }
   // restore the first block
   xx = full_rhs->begin ();
   for (auto kk = 0u; kk < blocks_size; kk++, xx++)
     *xx = *xx * rprec[kk];

  return (retval);
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
