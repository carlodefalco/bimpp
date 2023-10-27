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
    n = N - (size - 1)*N;
  }

  distributed_sparse_matrix A;
  A.set_ranges (n);

  for (int ir = 0; ir < size; ++ir) {
    if (ir == rank) {
      std::cout << "rank " << rank << std::endl;
      std::cout << "N " << N << std::endl;
      std::cout << "n " << n << std::endl;
      std::cout << "A.rows () " << A.rows () << std::endl;
      std::cout << "A.cols () " << A.cols () << std::endl;
      std::cout << "A.nnz () " << A.nnz () << std::endl;
    }
    MPI_Barrier (MPI_COMM_WORLD);
  }
  MPI_Finalize ();
  return (0);
}
