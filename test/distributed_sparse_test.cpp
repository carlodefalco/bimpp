/*
  Copyright (C) 2023 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <bim_sparse_distributed.h>
#include <bim_distributed_vector.h>

#include <mpi.h>
#include <fstream>

constexpr int N = 1000;
int main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  auto n = N / size;
  if (rank == size-1) {
    n = N - (size - 1)*n;
  } 
  
  distributed_sparse_matrix A;
  A.set_ranges (n);


  int start_row = A.range_start () < 5 ? A.range_start () : A.range_start () - 5;
  int end_row   = A.range_end () > N - 5 ? A.range_end () : A.range_end () + 5;

  for (int ii = start_row; ii < end_row; ++ii) {
    if (ii >= A.range_start () && ii < A.range_end ()) {
      A[ii][ii] = 4.;
      if (ii >= 5) {
	A[ii][ii-5] = -1.;
      }
      if (ii < N - 5) {
	A[ii][ii+5] = -1.;
      }
    } else {
      A[ii][ii] = 0.;
      if (ii >= 5) {
	A[ii][ii-5] = 0.;
      }
      if (ii < N - 5) {
	A[ii][ii+5] = 0.;
      }
    }
  }

  A.assemble ();

  distributed_vector x(n, MPI_COMM_WORLD);
  for (auto & ii : x.get_owned_data ()) {
    ii = 1.;
  }

  auto y = A*x;
  y. assemble ();
  
  for (int ir = 0; ir < size; ++ir) {
    if (ir == rank) {
      for (auto const & ii : y.get_owned_data ()) {
	std::cout << ii << std::endl;
      }
    
    }
    MPI_Barrier (MPI_COMM_WORLD);
  }
  MPI_Finalize ();
  return (0);
}
