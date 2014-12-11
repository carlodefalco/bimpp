
/*
	Problem:	du/dt-nabla(u)=g
						u=1-x^2-y^2-z^2 on border
						g=6

	Exact Solution:	u=1-x^2-y^2-z^2
*/

#include <stdio.h>
#include "lis_config.h"
#include "lis.h"

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <lis_operators.h>
#include <fstream>
#include <stdlib.h>

LIS_INT main(LIS_INT argc, char* argv[])
{

	int int_rank,rank;
	int int_size,size;
	
  LIS_INT     iter;
  double      time;
 
  lis_initialize(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD,&int_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	rank=int_rank;
	
	sparse_matrix       lhs,lhs_new,lhs_loc;
	std::vector<double> rhs1,rhs2,rhs_new,rhs_loc;
	
	std::vector<double> exactsolution_start, exactsolution,sol;

	std::vector<int> bnodes;
	std::vector<double> vnodes_start, vnodes;
	
	std::vector<double> uold;
	double dt=1;
	int T=10;

	Linear_Solver_Option option={"",0,""};

	std::ofstream fout_sol ("Lis_Solution_TimeTest.txt");

	if (rank == 0)
		{
			std::cout << "\n\n*****\nLis test: Time Dependent Problem\n*****\n";

		  std::cout << "read mesh" << std::endl;
		  mesh msh(std::string("mesh_in_cube.msh"));

		  std::cout << "compute mesh props" << std::endl;
		  msh.precompute_properties ();

		  std::cout << "assemble stiffness matrix. nnodes = " << msh.nnodes << std::endl;      

			bim3a_structure (msh, lhs);
		  std::vector<double> ecoeff (msh.nelements, 1.0); //isotropic diffusion coefficient
		  std::vector<double> v (msh.nnodes, 0.0);
			std::vector<double> ncoeff(msh.nnodes, 1/dt);
			std::vector<double> nodecoeff1(msh.nnodes,1/dt);
			std::vector<double> nodecoeff2(msh.nnodes,6.0);

			bim3a_advection_diffusion (msh, ecoeff, v, lhs);
		  bim3a_reaction (msh, ecoeff, ncoeff, lhs);

			bim3a_rhs (msh, ecoeff, nodecoeff1, rhs1);
			bim3a_rhs (msh, ecoeff, nodecoeff2, rhs2);
		
		  uold=std::vector<double>(msh.nnodes,0.0);
			for (int i=0;i<msh.nnodes;++i)
				{
					uold[i]=1-msh.p(0,i)*msh.p(0,i)-msh.p(1,i)*msh.p(1,i)-msh.p(2,i)*msh.p(2,i);								
				}      			

			std::vector<int> sidelist;
			sidelist.push_back(1);
			sidelist.push_back(2);
			sidelist.push_back(3);
			sidelist.push_back(4);
			sidelist.push_back(5);
			sidelist.push_back(6);

			bim3a_boundary_nodes(msh,sidelist,bnodes);
			vnodes_start.resize(bnodes.size());
			vnodes.resize(bnodes.size());
			for (int i=0;i<vnodes.size();++i)
				{	
					vnodes_start[i]=1-msh.p(0,bnodes[i])*msh.p(0,bnodes[i])-msh.p(1,bnodes[i])*msh.p(1,bnodes[i])
													-msh.p(2,bnodes[i])*msh.p(2,bnodes[i]);
				}
			exactsolution_start.resize (msh.nnodes);
			exactsolution.resize (msh.nnodes);
			for(int i=0;i<exactsolution_start.size();++i)
				{
					exactsolution_start[i]=1-msh.p(0,i)*msh.p(0,i)-msh.p(1,i)*msh.p(1,i)-msh.p(2,i)*msh.p(2,i);
				}
		
			std::cout << "\nResult of Time Dependent Test \nwill be written in output_time.txt\n";
		}

  for (int t=1;t<=T;++t)
		{
			if (rank == 0)
		    {  
					std::cout << "\nTime Iteration: "<<t<<std::endl;		
					rhs_new.resize(rhs1.size());
					lhs_new=lhs;

					for (int i=0;i<rhs_new.size();++i)
						{
							rhs_new[i]=rhs2[i]+rhs1[i]*uold[i];
						}
					for (int i=0;i<vnodes.size();++i)
						{	
							vnodes[i]=vnodes_start[i];
						}
							
					bim3a_dirichletBC(lhs_new,rhs_new,bnodes,vnodes);
								
					for(int i=0;i<exactsolution.size();++i)
						{
							exactsolution[i]=exactsolution_start[i];						
						}
					}
			if(size > 1)
				{
					lis_matrix_parallelization(lhs_new,lhs_loc);
					lis_vector_parallelization(rhs_new,rhs_loc);
						
					lis_solve_system(lhs_loc,rhs_loc,sol,iter,time,lhs_new.size(),option);
				}
			else
				lis_solve_system(lhs_new,rhs_new,sol,iter,time,lhs_new.size(),option); 	
			if (rank == 0)
				{
					std::cout<<"Number of iterations = "<<iter<<std::endl;
					std::cout<<"Elapsed time = "<< time<<std::endl;
					
				  fout_sol << std::endl;
					fout_sol << "Time Iteration: "<<t<<std::endl;
					double norm=0;
					//double normexact=0;
				  for (int i = 0; i < sol.size (); ++i)
						{
					
					  	fout_sol << sol[i] << "  " << exactsolution[i]<< std::endl;
							norm+=(exactsolution[i]-sol[i])*(exactsolution[i]-sol[i]);
							//normexact+=(exactsolution[k])*(exactsolution[k]);
							uold[i]=sol[i];
						}
				  
					std::cout<<"Error: "<<norm<<std::endl<<std::endl;
					assert(norm < 10e-2);
				}
		}
	fout_sol.close ();
  
  lis_finalize();
 
  return 0;
}

