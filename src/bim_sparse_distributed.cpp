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

  std::vector<int> temp1;
  std::vector<double> temp2;

  non_local.col_ind.swap (temp1);
  non_local.a.swap (temp2);
  
  this->set_properties ();

  non_local.prc_ptr.assign (this->mpisize + 1, 0);
  sparse_matrix::col_iterator jj;

  for (size_t ii = 0; ii < this->mpisize; ++ii)
    {
      non_local.prc_ptr[ii+1] = non_local.prc_ptr[ii];
      if (ii != mpirank)
        for (size_t kk = ranges[ii]; kk < ranges[ii+1]; ++kk)
          {
            // std::cout <<  non_local.prc_ptr.size () << " " << ii + 1 << std::endl;
            non_local.prc_ptr[ii+1] += (*this)[kk].size ();
            non_local.col_ind.reserve ((*this)[kk].size ());
            non_local.row_ind.reserve ((*this)[kk].size ());
            non_local.a.reserve ((*this)[kk].size ());
            for (jj  = (*this)[kk].begin ();
                 jj != (*this)[kk].end (); ++jj)
              {
                non_local.row_ind.push_back (kk);
                non_local.col_ind.push_back (this->col_idx (jj));
                non_local.a.push_back (this->col_val (jj));
              }
          }
    }
}


void
distributed_sparse_matrix::non_local_csr_update ()
{

  int idx = 0;
  for (int i = 1; i < non_local.row_ptr.size (); ++i)
    for (int j = 0; j < non_local.row_ptr[i] - non_local.row_ptr[i-1]; ++j)
      {
	non_local.a[idx] = (*this)[i-1][non_local.col_ind[idx]];
	idx ++;
      }
}
void
distributed_sparse_matrix::remap ()
{
  non_local_csr ();

  if (mpirank == 0) {
  std::cout << "rank " << mpirank << " non_local.prc_ptr : " << std::endl;
  for (auto ii : non_local.prc_ptr)
    std::cout << ii << " ";
  std::cout << std::endl;

  std::cout << "non_local.row_ind : " << std::endl;
  for (auto ii : non_local.row_ind)
    std::cout << ii << " ";
  std::cout << std::endl;

  std::cout << "non_local.col_ind : " << std::endl;
  for (auto ii : non_local.col_ind)
    std::cout << ii << " ";
  std::cout << std::endl;

  std::cout << "non_local.a : " << std::endl;
  for (auto ii : non_local.a)
    std::cout << ii << " ";
  std::cout << std::endl;
  }

  /// Distribute buffer sizes
  rank_nnz.assign (mpisize, 0);
  for (int ii = 0; ii < mpisize; ++ii)
    rank_nnz[ii] = non_local.prc_ptr[ii+1] - non_local.prc_ptr[ii];
  MPI_Alltoall (MPI_IN_PLACE, 1, MPI_INT, &(rank_nnz[0]), 1, MPI_INT, comm);

  /// Allocate buffers
  for (int ii = 0; ii < mpisize; ++ii)
    if ((rank_nnz[ii] > 0) && (ii != mpirank))
      {
        row_buffers[ii].resize (rank_nnz[ii]);
        col_buffers[ii].resize (rank_nnz[ii]);
        val_buffers[ii].resize (rank_nnz[ii]);
      }


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
      int rank_nnz_snd_ii = non_local.prc_ptr[ii+1] - non_local.prc_ptr[ii];
      if (rank_nnz_snd_ii > 0) // we must send something to rank ii
        {
          int send_tag = mpirank + mpisize * ii;
          reqs.resize (reqs.size () + 1);
          MPI_Isend (&(non_local.row_ind[non_local.prc_ptr[ii]]),
                     rank_nnz_snd_ii, MPI_INT, ii, send_tag, comm,
                     &(reqs.back ()));
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
      int rank_nnz_snd_ii = non_local.prc_ptr[ii+1] - non_local.prc_ptr[ii];
      if (rank_nnz_snd_ii > 0) // we must send something to rank ii
        {
          int send_tag = mpirank + mpisize * ii;
          reqs.resize (reqs.size () + 1);
          MPI_Isend (&(non_local.col_ind[non_local.prc_ptr[ii]]),
                     rank_nnz_snd_ii, MPI_INT, ii, send_tag, comm,
                     &(reqs.back ()));
        }     
    }
  MPI_Waitall (reqs.size (), &(reqs[0]), MPI_STATUSES_IGNORE);
  reqs.clear ();

  mapped = true;
}

void
distributed_sparse_matrix::assemble ()
{

  non_local_csr_update ()
  if (! mapped)
    remap ();

  /// 3) communicate values
  std::vector<MPI_Request> reqs;
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
      int rank_nnz_snd_ii = non_local.prc_ptr[ii+1] > non_local.prc_ptr[ii];      
      if (rank_nnz_snd_ii > 0) // we must send something to rank ii
        {
          int send_tag = mpirank + mpisize * ii;
          reqs.resize (reqs.size () + 1);
          MPI_Isend (&(non_local.a[non_local.prc_ptr[ii]]),
                     rank_nnz_snd_ii, MPI_DOUBLE, ii, send_tag, comm,
                     &(reqs.back ()));
        }     
    }
  MPI_Waitall (reqs.size (), &(reqs[0]), MPI_STATUSES_IGNORE);
  reqs.clear ();

  /// 4) insert communicated values into sparse_matrix
  for (int ii = 0; ii < mpisize; ++ii) // loop over ranks
    if (ii != mpirank)
      for (int kk = 0; kk < rank_nnz[ii]; ++kk)      
          (*this)[row_buffers[ii][kk]][col_buffers[ii][kk]]
            += val_buffers[ii][kk];


  /// 5) zero out communicated values  
  for (int iprc = 0; iprc < mpisize; ++iprc)
    if (iprc != mpirank)
      {
        for (int ii = non_local.prc_ptr[iprc];
             ii < non_local.prc_ptr[iprc+1]; ++ii)
          // the following will throw if we try to access
          // an element which does not exist yet!
          (*this)[non_local.row_ind[ii]].at (non_local.col_ind[ii]) = 0.0;
      }
}
