/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/
/*
  Problem:	-nabla(u)=g
  u=1-x^2-y^2-z^2 on border
  g=6

  Exact Solution:	u=1-x^2-y^2-z^2
*/

#include <stdio.h>
#include <lis_config.h>
#include <lis.h>
#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <lis_operators.h>
#include <fstream>
#include <stdlib.h>
#include <bim_config.h>

LIS_INT main(LIS_INT argc, char* argv[])
{
  int rank, size;

  LIS_INT     iter;
  double      time;

  lis_initialize(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  mesh msh;

  sparse_matrix       lhs, lhs_loc;
  std::vector<double> rhs, rhs_loc;

  std::vector<double> exactsolution,sol;
  linear_solver_option option = {"", "", ""};

  if (rank == 0)
    {  
      std::cout << "\n\n*****\nLis Test: Stationary Problem\n*****\n";

      std::cout << "read mesh" << std::endl;
      msh.read (data_dir + std::string("mesh_in_cube.msh"));

      std::cout << "compute mesh props" << std::endl;
      msh.precompute_properties ();

      bim3a_structure (msh, lhs);
      std::vector<double> ecoeff (msh.nelements, 1.0); //isotropic diffusion coefficient
      std::vector<double> v (msh.nnodes, 0.0);
      std::vector<double> ncoeff (msh.nnodes, 6.0);

      bim3a_advection_diffusion (msh, ecoeff, v, lhs);

      bim3a_rhs (msh, ecoeff, ncoeff, rhs);

      std::vector<int> sidelist;
      sidelist.push_back(1);
      sidelist.push_back(2);
      sidelist.push_back(3);
      sidelist.push_back(4);
      sidelist.push_back(5);
      sidelist.push_back(6);

      std::vector<int> bnodes;
      std::vector<double> vnodes;

      bim3a_boundary_nodes (msh, sidelist, bnodes);
      vnodes.resize (bnodes.size ());

      for(int i = 0; i < vnodes.size (); ++i)
	{
	  vnodes[i]=1-msh.p(0,bnodes[i])*msh.p(0,bnodes[i])-msh.p(1,bnodes[i])*msh.p(1,bnodes[i])-msh.p(2,bnodes[i])*msh.p(2,bnodes[i]);
	}

      bim3a_dirichletBC (lhs, rhs, bnodes, vnodes);

      exactsolution.resize (msh.nnodes);

      for(int i = 0; i < exactsolution.size (); ++i)
	{
	  exactsolution[i]=1-msh.p(0,i)*msh.p(0,i)-msh.p(1,i)*msh.p(1,i)-msh.p(2,i)*msh.p(2,i);
	}
    }
  if(size > 1)
    {
      lis_matrix_parallelization(lhs,lhs_loc);
      lis_vector_parallelization(rhs,rhs_loc);

      lis_solve_system(lhs_loc,rhs_loc,sol,iter,time,lhs.size(),option);
    }
  else
    lis_solve_system(lhs,rhs,sol,iter,time,lhs.size(),option);
  if (rank == 0)
    {
      std::cout<<"\nNumber of iterations = "<<iter<<std::endl;
      std::cout<<"Elapsed time = "<< time<<std::endl;
      std::cout << "\nResult of Stationary Test \nwill be written in output.txt\n";
      std::ofstream fout ("Lis_Solution.txt");
      fout << std::endl;

      double norm = 0;
      std::vector<double> delta (sol.size ());

      for (int k = 0; k < sol.size (); ++k)
        {
          fout << sol[k] << "  " << exactsolution[k]<< std::endl;
          delta[k]=exactsolution[k]-sol[k];
        }
      fout.close ();

      bim3a_norm (msh, delta, norm, Inf);
      std::cout << "Error: " << norm << std::endl;

      if (norm > 10e-10)
        {
          std::cerr << "The error is bigger than tolerance" << std::endl;
          exit(-1);
        }
    }

  lis_finalize();

  return 0;
}

