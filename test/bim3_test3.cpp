/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>
#include <bim_config.h>

int main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  if (rank == 0)
    {
      std::cout << "\n\n*****\ntest 3: Advection-Diffusion"
                << "  with diagonal diffusion tensor\n*****\n";

      std::cout << "read mesh" << std::endl;
      mesh msh (data_dir + std::string ("mesh_in.msh"));

      std::cout << "export mesh" << std::endl;
      msh.write (std::string ("mesh_out.m"));

      std::cout << "compute mesh props" << std::endl;
      msh.precompute_properties ();

      std::cout << "assemble stiffness matrix. nnodes = "
                << msh.nnodes << std::endl;
      sparse_matrix M;

      bim3a_structure (msh, M);
      std::vector<double> ecoeff (msh.nelements, 1.0);
      std::vector<double> v (msh.nnodes, 0.0);
      std::vector<double> b (msh.nnodes, 0.0);
      std::vector<double> ncoeff (msh.nnodes, 1.0);

      //diagonal diffusion tensor
      std::vector<double> dcoeff (msh.nelements * 3, 1.0);
      std::cout << "make advection-diffusion "
                << "with anisotropic diffusion" << std::endl;
      bim3a_advection_diffusion_anisotropic (msh, dcoeff, v, M);

      bim3a_reaction (msh, ecoeff, ncoeff, M);
      bim3a_rhs (msh, ecoeff, ncoeff, b);

      std::cout << "export stiffness matrix. nnodes = "
                << msh.nnodes << std::endl;
      std::ofstream fout ("SG3.m");
      fout << M;
      fout.close ();
    }

  MPI_Finalize ();
  return (0);
}
