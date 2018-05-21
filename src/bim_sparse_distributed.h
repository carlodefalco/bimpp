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
  
  void
    non_local_csr_update ();

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
  int nnz_local;
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
    update_assemble ();
  
  void
    remap ();


  void csr (std::vector<double> &a, std::vector<int> &col_ind,
	    std::vector<int> &row_ptr,
	    int base)
  {
    //this->set_properties ();

    a.resize (nnz_local); col_ind.resize (nnz_local);
    row_ptr.resize (ie - is + 1);
    int idx = 0;
    int idr = 0;

    typename sparse_matrix_template<double>::col_iterator jj;
    for (size_t ii = is; ii < ie; ++ii)
      {
	row_ptr[idr] = idx + base;

	if ((*this)[ii].size () > 0)
	  {
	    for (jj  = (*this)[ii].begin (); jj != (*this)[ii].end (); ++jj)
	      {
		col_ind[idx] = this->col_idx (jj) + base;
		a[idx] = this->col_val (jj);
		idx++;
	      }
	  }
	idr++;

      }

    std::fill (row_ptr.begin () + idr, row_ptr.end (), idx + base);
 
  }



  void csr_update (std::vector<double> &a,
		   const std::vector<int> &col_ind,
		   const std::vector<int> &row_ptr,
		   int base)
  {
    size_t ni = row_ptr.size ();
    size_t nj = col_ind.size ();
    a.clear ();
    a.reserve (nj);

    //std::cout << " ni = " << ni << std::endl;
    for (size_t in = 0; in < ni - 1; ++in)
      for (size_t jn = row_ptr[in] - base; jn < row_ptr[in+1] - base; ++jn)
	a.push_back (col_val (((*this)[in + is]).find (col_ind[jn] - base)));

  }

  
};

#endif
