/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#ifndef HAVE_BIM_SPARSE_DISTRIBUTED_H
#define HAVE_BIM_SPARSE_DISTRIBUTED_H 1

#include <mpi.h>
#include <bim_sparse.h>

class
distributed_sparse_matrix
  : public sparse_matrix
{

private :

  void
  non_local_csr ( );
  

  size_t is, ie;
  MPI_Comm comm;
  int mpirank, mpisize;

  struct
  non_local_t
  {
    std::vector<int> prc_ptr, row_ind, col_ind;
    std::vector<double> a;
  } non_local;

  std::map<int, std::vector<int>> row_buffers;
  std::map<int, std::vector<int>> col_buffers;
  std::map<int, std::vector<double>> val_buffers;

  std::vector<int> ranges;
  std::vector<int> rank_nnz;

  bool mapped;
  
public :

  void
  set_ranges (size_t is_, size_t ie_, MPI_Comm comm_ = MPI_COMM_WORLD);

  distributed_sparse_matrix (size_t is_, size_t ie_, MPI_Comm comm_ = MPI_COMM_WORLD)
    : mapped (false)
    { set_ranges (is_, ie_, comm_); }

  distributed_sparse_matrix (MPI_Comm comm_ = MPI_COMM_WORLD)
    : comm (comm_), mapped (false)
    { }

  void
  assemble ();

  void
  remap ();
  
};

#endif
