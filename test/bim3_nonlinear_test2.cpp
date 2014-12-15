
/*
	Problem:	-div(|grad(u)|^(p-2) grad(u))=f
						u=1/q*(0.5^q-[(x-0.5)^2+(y-0.5)^2+(z-0.5)^2)]^(q/2) on border
						f=3
						

	Exact Solution:	u=1/q*(0.5^q-[(x-0.5)^2+(y-0.5)^2+(z-0.5)^2)]^(q/2)
*/

#include <stdio.h>
#include "lis_config.h"
#include "lis.h"

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <lis_operators.h>
#include <nonlinear_solvers.h>
#include <nonlinear_operators.h>
#include <nonlinear_examples.h>
#include <fstream>
#include <stdlib.h>
#include <cmath>

LIS_INT main(LIS_INT argc, char* argv[])
{

	int rank;
	int size;
	
	lis_initialize(&argc, &argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
		
	std::vector<double> exactsolution, sol;

	std::vector<int> bnodes;
	std::vector<double> vnodes;

	std::vector<double> uold;
	
	std::vector<double> fcoeff;
	
	plaplacian plap;
	
	backtracking_inexact_newton_option option={100,10e-10, 10e-10, 0.5,10e-4,0.1,0.5,0,L2};
	
	stream_option stream={2,"Lis_Solution_NonLinear.txt"};
	linear_solver_option lis_option={"-tol 0.5","","-conv_cond 1 -i cg"};
	ForcingType1 forcing={0.9};

	inexact_newton_status status;

	double p=3;
	double q=p/(p-1);
	
	double lambda=1.0;
	double mu=0.5;
	
	if(rank==0)
		{
			std::cout << "\n\n*****\nLis test: Non Linear Problem\n*****\n";
						
			plap.read_mesh("mesh_in_cube.msh");
			int nnodes=plap.msh.nnodes;
	
			exactsolution.resize (nnodes);
			uold.resize(nnodes);
			for(int i=0;i<nnodes;++i)
				{
					exactsolution[i]=(1.0/q)*(pow(0.5,q)-pow((plap.msh.p(0,i)-0.5)*(plap.msh.p(0,i)-0.5)+(plap.msh.p(1,i)-0.5)*(plap.msh.p(1,i)-0.5)+(plap.msh.p(2,i)-0.5)*(plap.msh.p(2,i)-0.5),q/2.0));
					uold[i]=exactsolution[i]*(1+lambda*(plap.msh.p(0,i)-mu)*(plap.msh.p(1,i)-mu)*(plap.msh.p(2,i)-mu));		
				}	
			fcoeff.resize(nnodes,3.0);
	
			std::vector<int> sidelist; 
			sidelist.push_back(1);
			sidelist.push_back(2);
			sidelist.push_back(3);
			sidelist.push_back(4);
			sidelist.push_back(5);
			sidelist.push_back(6);

			bim3a_boundary_nodes(plap.msh,sidelist,bnodes);
			vnodes.resize(bnodes.size());
			for (int i=0;i<vnodes.size();++i)
				{	
					uold[bnodes[i]]=(1.0/q)*(pow(0.5,q)-pow((plap.msh.p(0,bnodes[i])-0.5)*(plap.msh.p(0,bnodes[i])-0.5)+(plap.msh.p(1,bnodes[i])-0.5)*(plap.msh.p(1,bnodes[i])-0.5)+(plap.msh.p(2,bnodes[i])-0.5)*(plap.msh.p(2,bnodes[i])-0.5),q/2));
					vnodes[i]=0.0;
				}
			
			plap.import(p,exactsolution,fcoeff,vnodes,bnodes);
			
		}			
	
	status=backtracking_inexact_newton<plaplacian,ForcingType1,LIS>(plap,uold,sol,forcing,lis_option,option,stream);
	
	if(rank == 0)
		{
			if(status.converged)
				{
					std::vector<double> delta_exact(sol.size());
					double deltaNorm;
					for(int i=0;i<sol.size();++i)
						delta_exact[i]=sol[i]-exactsolution[i];
					bim3a_norm(plap.msh,delta_exact,deltaNorm,L2);
					std::cout<<"Converged!"<<std::endl;
					std::cout<<"Total Newton's Iterations: "<<status.iteration<<std::endl;
					std::cout<<"Residual Norm: "<<status.residual<<std::endl;
					std::cout<<"Error: "<<deltaNorm<<std::endl;
				}
			else
				std::cerr<<"Not Converged =("<<std::endl;
		}

  lis_finalize();
 
  return 0;
}

