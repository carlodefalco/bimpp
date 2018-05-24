#include <cassert>
#include <iostream>
#include <vector>
#include <mpi.h>


//!    The data transfer pattern for assembly/synchronization 
//!    is shown in the graph below, I suspect there might be 
//!    some unnecessary redundancy that could be eliminated.
//!
//!        non_local                           owned
//!         ^    |                            |     ^
//!         |    |                            |     |
//! 6) copy |    | 1) copy             3) add |     | 4) copy
//!         |    |                            |     |
//!         |    |                            |     |
//!         |    V         5) send/recv       V     |
//!                <-------------------------
//!         ghosts                            mirrors
//!                ------------------------->
//!                        2) send/recv 


class
distributed_vector
{

private:

  MPI_Comm comm;

  int mpirank;
  int mpisize;
  
  int owned_count;
  int is, ie;

  bool mapped;
  
  /// vector entries owned by the current node
  std::vector<double> owned_data;

  /// vector entries touched by the current node
  /// that belong somwhere else
  std::map<int, double> non_local_data;

  struct
  ghosts_t
  {
    std::vector<int> prc_ptr, row_ind, rank_nnz;
    std::vector<double> a;
  } ghosts;

  void
  ghost_csr ()
  {
    /*
    std::cout << "non_local_data : " << std::endl;
    for (auto ii : non_local_data)
      std::cout << ii.first << " " << ii.second << " ";
    std::cout << std::endl;
    */
    ghosts.prc_ptr.assign (this->mpisize + 1, 0);
    ghosts.rank_nnz.assign (this->mpisize, 0);
    
    ghosts.row_ind.reserve (non_local_data.size ());
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

  
  struct
  mirrors_t
  {
    std::vector<int> prc_ptr, row_ind, rank_nnz;
    std::vector<double> a;
  } mirrors;

  /// the entries that are owned by rank i
  /// are numbered between ranges[i]
  /// and ranges[i+1]
  std::vector<int> ranges;

  /// check whether the idx-th global
  /// entry is owned by the current rank
  inline bool
  is_owned (int idx)
  { return idx >= is && idx < ie; }

  /// check whether the idx-th global
  /// entry is owned by the given rank
  inline bool
  is_owned (int idx, int irank)
  { return idx >= ranges[irank] && idx < ranges[irank+1]; }

  /// return the rank that owns the idx-th entry
  inline int
  owner (int idx)
  {
    auto ir = std::lower_bound (ranges.begin (), ranges.end (), idx);
    return int ((ir - ranges.begin ()) - 1);
  }
  
public:

  distributed_vector (int owned_count_,
                      MPI_Comm comm_ = MPI_COMM_WORLD)    
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


  distributed_vector (int is_, int ie_,
                      MPI_Comm comm_ = MPI_COMM_WORLD)
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
  operator() (int idx)
  {

    assert ((idx >= 0) && (idx < ranges.back ()));
    
    return is_owned (idx)
      ? owned_data[idx-is]
      : non_local_data[idx];
    
  }
  
  void
  remap ()
  {

    /// Step 1 : Copy non_local_data into ghosts
    ghost_csr ();

    /*
    std::cout << "ghosts.prc_ptr : " << std::endl;
    for (auto ii : ghosts.prc_ptr)
      std::cout << ii << " ";
    std::cout << std::endl;

    std::cout << "ghosts.row_ind : " << std::endl;
    for (auto ii : ghosts.row_ind)
      std::cout << ii << " ";
    std::cout << std::endl;

    std::cout << "ghosts.a : " << std::endl;
    for (auto ii : ghosts.a)
      std::cout << ii << " ";
    std::cout << std::endl;
    */
    
    /// Step 2 : Send ghosts to mirrors
    
    /// 2.1 : Distribute buffer sizes
    mirrors.rank_nnz = ghosts.rank_nnz;

    /*
    std::cout << "rank " << mpirank << " mirrors.rank_nnz : " << " ";
    for (auto ii : mirrors.rank_nnz)
      std::cout << ii << " ";
    std::cout << std::endl;

    std::cout << "rank " << mpirank << " ghosts.rank_nnz : " << " ";
    for (auto ii : ghosts.rank_nnz)
      std::cout << ii << " ";
    std::cout << std::endl;
    */
    
    MPI_Alltoall (MPI_IN_PLACE, 1, MPI_INT,
                  &(mirrors.rank_nnz[0]), 1, MPI_INT, comm);

    /*
    std::cout << "rank " << mpirank << " mirrors.rank_nnz : ";
    for (auto ii : mirrors.rank_nnz)
      std::cout << ii << " ";
    std::cout << std::endl;
    */
    mirrors.prc_ptr.assign (mpisize + 1, 0);

    //std::cout << "rank " << mpirank << " mirrors.prc_ptr ";
    for (int ii = 0; ii < mpisize; ++ii)
      {
        mirrors.prc_ptr[ii+1] = mirrors.prc_ptr[ii] +
          mirrors.rank_nnz[ii];
        //std::cout << mirrors.prc_ptr[ii+1] << " ";
      }
    //std::cout << std::endl;
    
    mirrors.row_ind.resize (mirrors.prc_ptr.back ());
    mirrors.a.resize (mirrors.prc_ptr.back ());

    /// 2.2 : Send ghost indices and receive mirror indices    
    std::vector<MPI_Request> reqs;
    for (int ii = 0; ii < mpisize; ++ii)
      {

        if (ii == mpirank) continue; // No communication to self!
        
        if (mirrors.rank_nnz[ii] > 0) // we must receive something from rank ii
          {
            int recv_tag = ii   + mpisize * mpirank;
            reqs.resize (reqs.size () + 1);
            MPI_Irecv (&(mirrors.row_ind[mirrors.prc_ptr[ii]]),
                       mirrors.rank_nnz[ii], MPI_INT, ii, recv_tag,
                       comm, &(reqs.back ()));
          }

        if (ghosts.rank_nnz[ii] > 0) // we must send something to rank ii
          {
            int send_tag = mpirank + mpisize * ii;
            reqs.resize (reqs.size () + 1);
            MPI_Isend (&(ghosts.row_ind[ghosts.prc_ptr[ii]]),
                       ghosts.rank_nnz[ii], MPI_INT, ii, send_tag,
                       comm, &(reqs.back ()));
          }
        
      }
  
      MPI_Waitall (reqs.size (), &(reqs[0]), MPI_STATUSES_IGNORE);
      reqs.clear ();

      mapped = true;
  }

  void
  assemble ()
  {
    
    if (! mapped)
      remap ();

    /// 2.3 : Send ghosts data and receive into mirrors
    std::vector<MPI_Request> reqs;
    for (int ii = 0; ii < mpisize; ++ii)
      {

        if (ii == mpirank) continue; // No communication to self!
        
        if (mirrors.rank_nnz[ii] > 0) // we must receive something from rank ii
          {
            int recv_tag = ii   + mpisize * mpirank;
            reqs.resize (reqs.size () + 1);
            MPI_Irecv (&(mirrors.a[mirrors.prc_ptr[ii]]),
                       mirrors.rank_nnz[ii], MPI_DOUBLE, ii,
                       recv_tag, comm, &(reqs.back ()));
          }

        if (ghosts.rank_nnz[ii] > 0) // we must send something to rank ii
          {
            int send_tag = mpirank + mpisize * ii;
            reqs.resize (reqs.size () + 1);
            MPI_Isend (&(ghosts.a[ghosts.prc_ptr[ii]]),
                       ghosts.rank_nnz[ii], MPI_DOUBLE, ii,
                       send_tag, comm, &(reqs.back ()));
          }
        
      }
  
      MPI_Waitall (reqs.size (), &(reqs[0]), MPI_STATUSES_IGNORE);
      reqs.clear ();

      /// Step 3 : Add mirrors into owned_data
      for (int ii = 0; ii < mirrors.prc_ptr.back (); ++ii)
        (*this)(mirrors.row_ind[ii]) += mirrors.a[ii];

      /// Step 4 : Copy owned_data into mirrors 
      for (int ii = 0; ii < mirrors.prc_ptr.back (); ++ii)
        mirrors.a[ii] = (*this)(mirrors.row_ind[ii]);

      /// Step 5 : Send mirrors data and receive into ghosts
      for (int ii = 0; ii < mpisize; ++ii)
      {

        if (ii == mpirank) continue; // No communication to self!
        
        if (ghosts.rank_nnz[ii] > 0) // we must receive something from rank ii
          {
            int recv_tag = ii   + mpisize * mpirank;
            reqs.resize (reqs.size () + 1);
            MPI_Irecv (&(ghosts.a[ghosts.prc_ptr[ii]]),
                       ghosts.rank_nnz[ii], MPI_DOUBLE, ii,
                       recv_tag, comm, &(reqs.back ()));
          }

        if (mirrors.rank_nnz[ii] > 0) // we must send something to rank ii
          {
            int send_tag = mpirank + mpisize * ii;
            reqs.resize (reqs.size () + 1);
            MPI_Isend (&(mirrors.a[mirrors.prc_ptr[ii]]),
                       mirrors.rank_nnz[ii], MPI_DOUBLE, ii,
                       send_tag, comm, &(reqs.back ()));
          }
        
      }
  
      MPI_Waitall (reqs.size (), &(reqs[0]), MPI_STATUSES_IGNORE);
      reqs.clear ();

      /// Step 6 : Copy ghosts data into non_local_data
      for (int ii = 0; ii < ghosts.prc_ptr.back (); ++ii)
        (*this)(ghosts.row_ind[ii]) = ghosts.a[ii];
  }


  
  friend std::ostream&
  operator<< (std::ostream &, distributed_vector&);
  
};

std::ostream&
operator<< (std::ostream &stream,
            distributed_vector& dv)
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
  
int
main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);
  
  int rank = 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  int size = 0;
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  distributed_vector dv (10);

  if (rank == 0)
    {
      for (int ii = 0; ii < 13; ++ii)
        dv(ii) = 1;
    }

  if (rank == 1)
    {
      for (int ii = 8; ii < 20; ++ii)
        dv(ii) = -2;
    }

  dv.assemble ();
  MPI_Barrier (MPI_COMM_WORLD);

  for (int irank = 0; irank < size; ++irank)
    {
      if (rank == irank)
        {
          std::cout << "rank : " << rank << std::endl;
          std::cout << dv;
        }
      MPI_Barrier (MPI_COMM_WORLD);
    }

  
  MPI_Finalize ();
  return 0;
}
