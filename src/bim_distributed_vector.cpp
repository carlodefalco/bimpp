/*
  Copyright (C) 2018,2019 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/


#include <bim_distributed_vector.h>
#include <mpi_alltoallv2.h>

#include <cassert>
#include <iostream>

binary_operator replace_op =
  [] (const double & x, const double & y)
{
  return x;
};

binary_operator max_op =
  [] (const double & x, const double & y)
{
  return std::max (x, y);
};

void
distributed_vector::ghost_csr ()
{

  ghosts.prc_ptr.assign (this->mpisize + 1, 0);
  ghosts.rank_nnz.assign (this->mpisize, 0);
  ghosts.row_ind.clear();
  ghosts.row_ind.reserve (non_local_data.size ());
  ghosts.a.clear();
  ghosts.a.reserve (non_local_data.size ());

    
  int irank = 0;
  for (auto ii = non_local_data.begin ();
       (ii != non_local_data.end ())
         && irank < mpisize;
       ++ii)
    {
        
      while (! is_owned (ii->first, irank)) ++irank;

      (ghosts.rank_nnz[irank])++;
      ghosts.row_ind.push_back (ii->first);
      ghosts.a.push_back (ii->second);
        
    }    
    
  for (int ii = 0; ii < mpisize; ++ii)
    ghosts.prc_ptr[ii+1] = ghosts.prc_ptr[ii] +
      ghosts.rank_nnz[ii];
    
}

void
distributed_vector::ghost_csr_update ()
{
  int ind = 0;
  for (auto ii = non_local_data.begin ();
       ii != non_local_data.end ();
       ++ii)
    ghosts.a[ind++] = (ii->second);
}

distributed_vector::distributed_vector (int owned_count_,
                                        MPI_Comm comm_)    
  : owned_count (owned_count_), comm (comm_), mapped (false) 
{
    
  MPI_Comm_size (comm, &mpisize);
  MPI_Comm_rank (comm, &mpirank);

  owned_data.assign (owned_count, .0);

  /// Gather ranges
  ranges.assign (mpisize + 1, 0);
  MPI_Allgather (&owned_count, 1, MPI_INT, &(ranges[1]),
                 1, MPI_INT, comm);    
  for (int irank = 0; irank < mpisize; ++irank)
    ranges[irank+1] += ranges[irank];

  is = ranges[mpirank];
  ie = ranges[mpirank+1];
    
}


distributed_vector::distributed_vector (int is_, int ie_,
                                        MPI_Comm comm_)
  :  is (is_), ie (ie_), comm (comm_), mapped (false) 
{
    
  MPI_Comm_size (comm, &mpisize);
  MPI_Comm_rank (comm, &mpirank);

  owned_data.assign (ie - is, .0);

  /// Gather ranges
  ranges.assign (mpisize + 1, 0);
  MPI_Allgather (&ie, 1, MPI_INT, &(ranges[1]),
                 1, MPI_INT, comm);    

  owned_count = ie - is;
    
}

double&
distributed_vector::operator() (int idx)
{
  assert ((idx >= 0) && (idx < this->ranges.back ()));
    
  return is_owned (idx)
    ? owned_data[idx-is]
    : non_local_data[idx];
    
}
  
const double&
distributed_vector::operator() (int idx) const
{
  assert ((idx >= 0) && (idx < this->ranges.back ()));
    
  return is_owned (idx)
    ? owned_data[idx-is]
    : non_local_data.at (idx);
}
  
void
distributed_vector::remap ()
{
  /// Step 1 : Copy non_local_data into ghosts
  ghost_csr ();
    
  /// Step 2 : Send ghosts to mirrors
    
  /// 2.1 : Distribute buffer sizes
  mirrors.rank_nnz = ghosts.rank_nnz;
    
  MPI_Alltoall (MPI_IN_PLACE, 1, MPI_INT,
                &(mirrors.rank_nnz[0]), 1, MPI_INT, comm);

  mirrors.prc_ptr.assign (mpisize + 1, 0);

  for (int ii = 0; ii < mpisize; ++ii)
      mirrors.prc_ptr[ii+1] = mirrors.prc_ptr[ii] +
        mirrors.rank_nnz[ii];
    
  mirrors.row_ind.resize (mirrors.prc_ptr.back ());
  mirrors.a.resize (mirrors.prc_ptr.back ());

  /// 2.2 : Send ghost indices and receive mirror indices    
  void * sendbuf[mpisize];
  int sendcnts[mpisize];
  void * recvbuf[mpisize];
  int recvcnts[mpisize];

  for (int ii = 0; ii < mpisize; ++ii)
    {
      recvbuf[ii]  = &(mirrors.row_ind[mirrors.prc_ptr[ii]]);
      recvcnts[ii] =   mirrors.rank_nnz[ii];
      sendbuf[ii]  = &(ghosts.row_ind[ghosts.prc_ptr[ii]]);
      sendcnts[ii] =   ghosts.rank_nnz[ii];
    }

  MPI_Alltoallv2 (sendbuf, sendcnts, MPI_INT,
                  recvbuf, recvcnts, MPI_INT, comm);

  mapped = true;
}

void
distributed_vector::assemble (const binary_operator & binary_op)
{
  if (! mapped)
    remap ();
  else
    ghost_csr_update ();
  
  /// 2.3 : Send ghosts data and receive into mirrors
  void * sendbuf[mpisize];
  int sendcnts[mpisize];
  void * recvbuf[mpisize];
  int recvcnts[mpisize];

  for (int ii = 0; ii < mpisize; ++ii)
    {
      recvbuf[ii]  = &(mirrors.a[mirrors.prc_ptr[ii]]);
      recvcnts[ii] =   mirrors.rank_nnz[ii];
      sendbuf[ii]  = &(ghosts.a[ghosts.prc_ptr[ii]]);
      sendcnts[ii] =   ghosts.rank_nnz[ii];
    }

  MPI_Alltoallv2 (sendbuf, sendcnts, MPI_DOUBLE,
                  recvbuf, recvcnts, MPI_DOUBLE, comm);

  /// Step 3 : Add mirrors into owned_data
  for (int ii = 0; ii < mirrors.prc_ptr.back (); ++ii)
    (*this)(mirrors.row_ind[ii]) =
      binary_op ((*this)(mirrors.row_ind[ii]),
                 mirrors.a[ii]);

  /// Step 4 : Copy owned_data into mirrors 
  for (int ii = 0; ii < mirrors.prc_ptr.back (); ++ii)
    mirrors.a[ii] = (*this)(mirrors.row_ind[ii]);

  /// Step 5 : Send mirrors data and receive into ghosts
  for (int ii = 0; ii < mpisize; ++ii)
    {
      recvbuf[ii]  = &(ghosts.a[ghosts.prc_ptr[ii]]);
      recvcnts[ii] =   ghosts.rank_nnz[ii];
      sendbuf[ii]  = &(mirrors.a[mirrors.prc_ptr[ii]]);
      sendcnts[ii] =   mirrors.rank_nnz[ii];
    }
  
  MPI_Alltoallv2 (sendbuf, sendcnts, MPI_DOUBLE,
                  recvbuf, recvcnts, MPI_DOUBLE, comm);
  
  /// Step 6 : Copy ghosts data into non_local_data
  for (int ii = 0; ii < ghosts.prc_ptr.back (); ++ii)
    (*this)(ghosts.row_ind[ii]) = ghosts.a[ii];
}

void
distributed_vector::clear_non_local ()
{
  non_local_data.clear ();
  
  ghosts.prc_ptr.clear ();
  ghosts.row_ind.clear ();
  ghosts.rank_nnz.clear ();
  ghosts.a.clear ();
  
  mapped = false;
}

std::ostream&
operator<< (std::ostream &stream,
            const distributed_vector& dv)
{

  auto ii = dv.non_local_data.begin ();
  stream << "non_local_data : " << std::endl;
  for (;
       (ii != dv.non_local_data.end ())
         && (ii->first < dv.is);
       ++ii)
    stream << "idx = " << ii->first << " val = "
           << ii->second << std::endl;

  stream << "owned_data : " << std::endl;
  int jj = dv.is;
  for (auto kk : dv.owned_data)
    stream << "idx = " << jj++ << " val = "
           << kk << std::endl;

  stream << "non_local_data : " << std::endl;
  for (;
       (ii != dv.non_local_data.end ());
       ++ii)
    stream << "idx = " << ii->first << " val = "
           << ii->second << std::endl;
  return stream;
}

distributed_vector
operator * (distributed_sparse_matrix& M, const distributed_vector& x)
{
  distributed_vector y (x.get_range_start (), x.get_range_end ());
  distributed_vector X (x.get_range_start (), x.get_range_end ());
  X.get_owned_data () = x.get_owned_data ();
  for (auto ir = M.range_start (); ir < M.range_end (); ++ir) {
    for (auto jc = M[ir].begin (); jc != M[ir].end (); ++jc) {
      X(M.col_idx (jc)) += 0.;
    }
  }
  
  X.assemble (replace_op);
  for (auto ir = M.range_start (); ir < M.range_end (); ++ir) {
    for (auto jc = M[ir].begin (); jc != M[ir].end (); ++jc) {
      y(ir) += M.col_val (jc) * X[M.col_idx (jc)];
    }
  }
  return y;
}


distributed_vector
operator * (sparse_matrix& M, const distributed_vector& x)
{
  distributed_vector y (x.get_range_start (), x.get_range_end ());
  sparse_matrix::col_iterator j;
  for (unsigned int i = 0; i < M.size (); ++i)
    if (M[i].size ())
      for (j = M[i].begin (); j != M[i].end (); ++j)
        if (M.col_val (j) != 0)
          y[i] += M.col_val (j) * x[M.col_idx (j)];

  return y;
}

