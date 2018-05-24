/*
  Copyright (C) 2018 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/


#include <bim_distributed_vector.h>

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
