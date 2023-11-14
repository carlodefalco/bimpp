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

  distributed_vector x(n, MPI_COMM_WORLD);
  for (auto & ii : x.get_owned_data ()) {
    ii = 1.;
  }
  for (auto irank = 0; irank < size-1; ++irank) {
    if (rank == irank) {
      x(x.get_range_end()+3) = 3;
    }
  }

  x. assemble ();

  for (auto irank = 0; irank < size-1; ++irank) {
    if (rank == irank) {
      x[x.get_range_end()] = 3;
    }
  }
  x.remap();
  x.assemble();

  for (int ir = 0; ir < size; ++ir) {
    if (ir == rank) {
      for (auto const & ii : x.get_owned_data ()) {
	      std::cout << ii << std::endl;
      }
    
    }
    MPI_Barrier (MPI_COMM_WORLD);
  }
  MPI_Finalize ();
  return (0);
}
