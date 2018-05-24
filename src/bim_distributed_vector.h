/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#ifndef HAVE_BIM_DISTRIBUTED_VECTOR_H
#define HAVE_BIM_DISTRIBUTED_VECTOR_H 1

#include <iostream>
#include <map>
#include <mpi.h>
#include <vector>

//!    Class for distributed memory vector.
//|
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

  /// Structure to hold data of ghost entries
  /// in a format amenable for send/receive
  struct
  ghosts_t
  {
    std::vector<int> prc_ptr, row_ind, rank_nnz;
    std::vector<double> a;
  } ghosts;

  /// Copy data from non_local_data to ghosts
  void
  ghost_csr ();

  /// Structure to hold data of ghost entries
  /// in a format amenable for send/receive
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
                      MPI_Comm comm_ = MPI_COMM_WORLD);


  distributed_vector (int is_, int ie_,
                      MPI_Comm comm_ = MPI_COMM_WORLD);

  double&
  operator() (int idx);

  double&
  operator[] (int idx)
  { return (*this)(idx); };

  void
  remap ();

  void
  assemble ();

  friend std::ostream&
  operator<< (std::ostream &, distributed_vector&);
  
};

#endif

