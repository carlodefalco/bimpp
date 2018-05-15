/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/


#include <bim_sparse_distributed.h>

void
distributed_sparse_matrix::set_ranges (size_t is_,
                                       size_t ie_,
                                       MPI_Comm comm_)
{
  is = is_; ie = ie_; comm = comm_;
  MPI_Comm_rank (comm, &mpirank);
  MPI_Comm_size (comm, &mpisize);

  /// Gather ranges
  ranges.assign (mpisize + 1, 0);
  MPI_Allgather (&ie, 1, MPI_INT, &(ranges[1]), 1, MPI_INT, comm);

}


void
distributed_sparse_matrix::non_local_csr ()
{
  this->set_properties ();
  non_local.a.reserve (this->nnz / this->mpisize);
  non_local.col_ind.reserve (this->nnz / this->mpisize);
  non_local.row_ptr.resize (this->rows () + 1);

  int idx = 0, idr = 0;

  sparse_matrix::col_iterator jj;

  for (size_t ii = 0; ii < this->size (); ++ii)
    {
      non_local.row_ptr[idr++] = idx;
      if (((*this)[ii].size () > 0) && ((ii < is) || (ii >= ie)))
        {
          for (jj  = (*this)[ii].begin ();
               jj != (*this)[ii].end (); ++jj)
            {
              non_local.col_ind.push_back (this->col_idx (jj));
              non_local.a.push_back (this->col_val (jj));
              idx++;
            }
        }
    }
  std::fill (non_local.row_ptr.begin () + idr,
             non_local.row_ptr.end (),
             non_local.col_ind.size ());
}

void
distributed_sparse_matrix::remap ()
{
  non_local_csr ();
    
  /// Distribute buffer sizes
  rank_nnz.assign (mpisize, 0);
  for (int ii = 0; ii < mpisize; ++ii)
    rank_nnz[ii] = non_local.row_ptr[ranges[ii+1]]
      - non_local.row_ptr[ranges[ii]];
  MPI_Alltoall (MPI_IN_PLACE, 1,
                MPI_INT, &(rank_nnz[0]),
                1, MPI_INT, comm);

 
  /// Allocate buffers
  for (int ii = 0; ii < mpisize; ++ii)
    if ((rank_nnz[ii] > 0) && (ii != mpirank))
      {
        row_buffers[ii].resize (ranges[mpirank+1]
                                - ranges[mpirank] + 1);
        col_buffers[ii].resize (rank_nnz[ii]);
        val_buffers[ii].resize (rank_nnz[ii]);
      }
  mapped = true;
}

void
distributed_sparse_matrix::assemble ()
{

  if (! mapped)
    remap ();
  
  /// Communicate overlap regions
  
  /// 1) communicate row_ptr
  std::vector<MPI_Request> reqs;
  for (int ii = 0; ii < mpisize; ++ii)
    {
      if (ii == mpirank) continue; // No communication to self!
      if (rank_nnz[ii] > 0) // we must receive something from rank ii
        {
          int recv_tag = ii   + mpisize * mpirank;
          reqs.resize (reqs.size () + 1);
          MPI_Irecv (&(row_buffers[ii][0]), row_buffers[ii].size (),
                     MPI_INT, ii, recv_tag, comm, &(reqs.back ()));
        }
      if (non_local.row_ptr[ranges[ii+1]]
          > non_local.row_ptr[ranges[ii]]) // we must send something to rank ii
        {
          int send_tag = mpirank + mpisize * ii;
          reqs.resize (reqs.size () + 1);
          MPI_Isend (&(non_local.row_ptr[ranges[ii]]),
                     ranges[ii+1] - ranges[ii] + 1,
                     MPI_INT, ii, send_tag, comm, &(reqs.back ()));
        }     
      }
  
  MPI_Waitall (reqs.size (), &(reqs[0]), MPI_STATUSES_IGNORE);
  reqs.clear ();
    
  /// 2) communicate col_ind
  for (int ii = 0; ii < mpisize; ++ii)
    {
      if (ii == mpirank) continue; // No communication to self!
      if (rank_nnz[ii] > 0) // we must receive something from rank ii
        {
          int recv_tag = ii   + mpisize * mpirank;
          reqs.resize (reqs.size () + 1);
          MPI_Irecv (&(col_buffers[ii][0]), col_buffers[ii].size (),
                     MPI_INT, ii, recv_tag, comm, &(reqs.back ()));
        }
      if (non_local.row_ptr[ranges[ii+1]]
          > non_local.row_ptr[ranges[ii]]) // we must send something to rank ii
        {
          int send_tag = mpirank + mpisize * ii;
          reqs.resize (reqs.size () + 1);
          MPI_Isend (&(non_local.col_ind
                       [non_local.row_ptr[ranges[ii]]]),
                     non_local.row_ptr[ranges[ii+1]]
                     - non_local.row_ptr[ranges[ii]],
                     MPI_INT, ii, send_tag, comm, &(reqs.back ()));
        }     
      }
  MPI_Waitall (reqs.size (), &(reqs[0]), MPI_STATUSES_IGNORE);
  reqs.clear ();
    
  /// 3) communicate values
  for (int ii = 0; ii < mpisize; ++ii)
    {
      if (ii == mpirank) continue; // No communication to self!
      if (rank_nnz[ii] > 0) // we must receive something from rank ii
        {
          int recv_tag = ii   + mpisize * mpirank;
          reqs.resize (reqs.size () + 1);
          MPI_Irecv (&(val_buffers[ii][0]), val_buffers[ii].size (),
                     MPI_DOUBLE, ii, recv_tag, comm, &(reqs.back ()));
        }
      if (non_local.row_ptr[ranges[ii+1]]
          > non_local.row_ptr[ranges[ii]]) // we must send something to rank ii
        {
          int send_tag = mpirank + mpisize * ii;
          reqs.resize (reqs.size () + 1);
          MPI_Isend (&(non_local.a[non_local.row_ptr[ranges[ii]]]),
                     non_local.row_ptr[ranges[ii+1]]
                     - non_local.row_ptr[ranges[ii]],
                     MPI_DOUBLE, ii, send_tag, comm, &(reqs.back ()));
        }     
      }
  MPI_Waitall (reqs.size (), &(reqs[0]), MPI_STATUSES_IGNORE);
  reqs.clear ();
  
  /// 4) insert communicated values into sparse_matrix
  for (int ii = 0; ii < mpisize; ++ii) // loop over ranks
    if ((rank_nnz[ii] > 0) && (ii != mpirank))
      for (int jj = 0;
           jj < ranges[mpirank+1] - ranges[mpirank];
           ++jj) // loop over rows of rank
        for (int kk = row_buffers[ii][jj] - row_buffers[ii][0];
             kk < row_buffers[ii][jj+1]   - row_buffers[ii][0];
             ++kk) // loop over columns of row
          (*this)[jj + is][col_buffers[ii][kk]]
            += val_buffers[ii][kk];

  /// 5) zero out communicated values  
  for (int irow = 0; irow < is; ++irow)
    if ((*this)[irow].size ())
      for (auto& jj : (*this)[irow])
        jj.second = 0.0;

  for (int irow = ie; irow < this->size (); ++irow)
    if ((*this)[irow].size ())
      for (auto& jj : (*this)[irow])
        jj.second = 0.0;
}
