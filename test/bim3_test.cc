/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/

#include <sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>
  
int main (void)
{

  MPI::Init ();
  int rank = MPI::COMM_WORLD.Get_rank ();
  int size = MPI::COMM_WORLD.Get_size ();
    
  if (rank == 0)
    {  
      std::cout << "\n\n*****\ntest 1\n*****\n";
      
      std::cout << "read mesh" << std::endl;
      mesh msh (std::string("mesh_in.msh"));

      std::cout << "export mesh" << std::endl;
      msh.write (std::string("mesh_out.m"));
      
      std::cout << "compute mesh props" << std::endl;
      msh.precompute_properties ();
      
      std::cout << "assemble stiffness matrix. nnodes = " << msh.nnodes << std::endl;      
      sparse_matrix M;
      
      bim3a_structure (msh, M);
      std::vector<double> ecoeff (msh.nelements, 1.0);
      std::vector<double> v (msh.nnodes, 0.0);
      std::vector<double> b (msh.nnodes, 0.0);
      std::vector<double> ncoeff (msh.nnodes, 1.0);
      
      bim3a_advection_diffusion (msh, ecoeff, v, M);
      bim3a_reaction (msh, ecoeff, ncoeff, M);
      bim3a_rhs (msh, ecoeff, ncoeff, b);
      
      std::cout << "export stiffness matrix. nnodes = " << msh.nnodes << std::endl;      
      std::ofstream fout ("SG.m");
      fout << M;
      fout.close ();
    }
  
  MPI::Finalize();
  return (0);
}
