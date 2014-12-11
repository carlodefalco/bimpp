#ifndef HAVE_NL_SOLVERS_H
#define HAVE_NL_SOLVERS_H 1
#include <lis.h>
#include <lis_operators.h>
#include <mumps_class.h>
#include <fstream>
#include <nonlinear_operators.h>
#include <string>

enum SolverType {MUMPS,LIS};

template <typename P, typename FT, SolverType ST>
Inexact_Newton_Status inexact_newton(P& problem,
										std::vector<double>& ustart,
										std::vector<double>& sol,
										FT& forcing,
										Linear_Solver_Option lis_option={"",0,"-conv_cond 1 -i cg","-tol 0.5"},
										Inexact_Newton_Option option={100,10e-10,10e-10,10e-10,Inf},
										const Stream_Option stream={2,"output.txt"})
{
	std::ofstream fout;
	
	double stepNorm=1;
	double resNorm=1;
	sparse_matrix lhs;
	std::vector<double> rhs,du,res; 
	sparse_matrix lhs_loc;
	std::vector<double> rhs_loc;
	std::vector<double> exactsolution;
	problem.get_solution(exactsolution);

	int rank,size;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	if(rank==0 && stream.verbosity>=1)
		{
			fout.open(stream.filename);
			std::cout << "Result of Non Linear Test \nwill be written in "<<stream.filename<<"\n\n";
		}
	if(rank==0)
		{
			problem(lhs,rhs,ustart);
			bim3a_norm(problem.msh,rhs,resNorm,option.type);
		}
	MPI_Bcast(&resNorm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if(resNorm<=option.minRes)
		{
			sol=ustart;
			return {0,resNorm,true};
		}

	unsigned int it(0);
	do
		{
			++it;
			if(rank==0 && stream.verbosity==2)		
				std::cout<<"Newton Iteration: "<<it<<std::endl;
						
			if(ST==LIS)
				{	
					LIS_INT iter;
					double time;
					if(size>1)
						{
							lis_matrix_parallelization(lhs,lhs_loc);
							lis_vector_parallelization(rhs,rhs_loc);
							lis_solve_system(lhs_loc,rhs_loc,du,iter,time,lhs.size(),lis_option);
						}
					else
						lis_solve_system(lhs,rhs,du,iter,time,lhs.size(),lis_option);
					if(rank==0 && stream.verbosity==2)
						{
							std::cout<<"Number of iterations = "<<iter<<std::endl;
							std::cout<<"Elapsed time = "<< time<<std::endl;
						}
				}
			else
				{	
					std::vector<double> xa;
					std::vector<int> ir,jc;
					if(rank == 0) 
						lhs.aij(xa,ir,jc,1);
						
					mumps mumps_solver;
  
					if (rank == 0)
						mumps_solver.set_lhs_structure (lhs.rows (), ir, jc);
		
					mumps_solver.analyze ();
		
					if (rank == 0)
						mumps_solver.set_lhs_data (xa);

					mumps_solver.factorize ();
		
					if (rank == 0)
						mumps_solver.set_rhs (rhs);

					mumps_solver.solve ();
				
					if (rank == 0)
						du=rhs;

					mumps_solver.cleanup ();
				}	
			if (rank == 0)
				{
					if(stream.verbosity==2)
						fout << "Iteration: "<<it<<std::endl;
						
					option.forcing=forcing(problem, ustart,du,option.forcing,option.type);
					
				  for (int i = 0; i < du.size (); ++i)
						ustart[i]=du[i]+ustart[i];
						
					bim3a_norm(problem.msh,du,stepNorm,option.type);
					
					problem(lhs,rhs,ustart);
					bim3a_norm(problem.msh,rhs,resNorm,option.type);
					if(stream.verbosity==2)
						{
							for (int i = 0; i < ustart.size (); ++i)
								fout << ustart[i]<<std::endl;
							
							//bim3a_print(ustart,fout);
							std::cout<<"Step Error: "<<stepNorm<<std::endl;
							std::cout<<"Residual Error: "<<resNorm<<std::endl<<std::endl;
						}
				}
			MPI_Bcast(&stepNorm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&resNorm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&option.forcing,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			sprintf(lis_option.tolerance,"-tol %f", option.forcing);
		}
		while(it<=option.maxIter && resNorm>option.minRes && stepNorm>option.tol);
		
		if(rank==0 && stream.verbosity==1)
			{
				fout<<"Solution: "<<std::endl;
				for (int i = 0; i < ustart.size (); ++i)
					fout << ustart[i]<<std::endl;
				
				//bim3a_print(ustart,fout);
			}
		if(rank==0 && stream.verbosity>=1)
			fout.close();

		bool converged(it<=option.maxIter);
		sol=ustart;
		return {it,resNorm,converged};
};

template <typename P, typename FT, SolverType ST>
Inexact_Newton_Status backtracking_inexact_newton(P& problem,
										std::vector<double>& ustart,
										std::vector<double>& sol,
										FT& forcing,
										Linear_Solver_Option lis_option={"",0,"-conv_cond 1 -i cg","-tol 0.5"},
										Backtracking_Inexact_Newton_Option option={100,10e-10,10e-10,10e-10,10e-4,0.1,0.5,0,Inf},
										const Stream_Option stream={2,"output.txt"})
{
	std::ofstream fout;
	
	double stepNorm=1;
	double resNorm=1;
	sparse_matrix lhs;
	std::vector<double> rhs,du,res; 
	sparse_matrix lhs_loc;
	std::vector<double> rhs_loc;
	std::vector<double> exactsolution;
	problem.get_solution(exactsolution);

	int rank,size;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	if(rank==0 && stream.verbosity>=1)
		{
			fout.open(stream.filename);
			std::cout << "Result of Non Linear Test \nwill be written in "<<stream.filename<<"\n\n";
		}
	if(rank==0)
		{
			problem(lhs,rhs,ustart);
			bim3a_norm(problem.msh,rhs,resNorm,option.type);
		}
	MPI_Bcast(&resNorm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if(resNorm<=option.minRes)
		{
			sol=ustart;
			return {0,resNorm,true};
		}

	unsigned int it(0);
	do
		{
			++it;
			if(rank==0 && stream.verbosity==2)		
				std::cout<<"Newton Iteration: "<<it<<std::endl;
						
			if(ST==LIS)
				{	
					LIS_INT iter;
					double time;
					if(size>1)
						{
							lis_matrix_parallelization(lhs,lhs_loc);
							lis_vector_parallelization(rhs,rhs_loc);
			
							lis_solve_system(lhs_loc,rhs_loc,du,iter,time,lhs.size(),lis_option);
						}
					else
						lis_solve_system(lhs,rhs,du,iter,time,lhs.size(),lis_option);
					if(rank==0 && stream.verbosity==2)
						{
							std::cout<<"Number of iterations = "<<iter<<std::endl;
							std::cout<<"Elapsed time = "<< time<<std::endl;
						}
				}
			else
				{	
					std::vector<double> xa;
					std::vector<int> ir,jc;
					if(rank == 0) 
						lhs.aij(xa,ir,jc,1);
						
					mumps mumps_solver;
  
					if (rank == 0)
						mumps_solver.set_lhs_structure (lhs.rows (), ir, jc);
		
					mumps_solver.analyze ();
		
					if (rank == 0)
						mumps_solver.set_lhs_data (xa);

					mumps_solver.factorize ();
		
					if (rank == 0)
						mumps_solver.set_rhs (rhs);

					mumps_solver.solve ();
				
					if (rank == 0)
						du=rhs;

					mumps_solver.cleanup ();
				}	


			if (rank == 0)
				{

					if(stream.verbosity==2)
						fout << "Iteration: "<<it<<std::endl;
					
					std::vector<double> unew(ustart.size());
				  for (int i = 0; i < du.size (); ++i)
						unew[i]=du[i]+ustart[i];
					
					std::vector<double> Fnew,Fold;
					double foldNorm=0;
					double fnewNorm=0;
					problem(Fold,ustart);
					problem(Fnew,unew);
					bim3a_norm(problem.msh,Fold,foldNorm,option.type);
					
					bim3a_norm(problem.msh,Fnew,fnewNorm,option.type);
					
					while(fnewNorm>(1-option.t*(1-option.forcing))*foldNorm)
						{
							bool flag=fnewNorm>(1-option.t*(1-option.forcing))*foldNorm;
							std::cout<<flag<<std::endl;
	
							std::vector<double> temp;
							double tempNorm=0;
							
							bim3a_matrix_vector_product(lhs,du,temp);
							for(int i=0;i<rhs.size();++i)
								tempNorm+=rhs[i]*temp[i];
		
							double a=fnewNorm*fnewNorm-foldNorm*foldNorm-2*tempNorm;
							double b=2*tempNorm;
							double c=foldNorm*foldNorm;
							
							option.theta_choice(a,b,c);
							
							for(int i=0;i<du.size();++i)
								{
									du[i]*=option.theta;
									unew[i]=du[i]+ustart[i];		
								}			
							option.forcing=1-option.theta*(1-option.forcing);
							
							problem(Fnew,unew);
							bim3a_norm(problem.msh,Fnew,fnewNorm,option.type);
						}
					option.forcing=forcing(problem,ustart,du,option.forcing,option.type);

					for (int i = 0; i < du.size (); ++i)
						ustart[i]=du[i]+ustart[i];	
					stepNorm=0;
					resNorm=0;
					bim3a_norm(problem.msh,du,stepNorm,option.type);
					
					problem(lhs,rhs,ustart);
					bim3a_norm(problem.msh,rhs,resNorm,option.type);
				
					if(stream.verbosity==2)
						{
							for (int i = 0; i < ustart.size (); ++i)
								fout << ustart[i]<<std::endl;
							
							//bim3a_print(ustart,fout);
							std::cout<<"Step Error: "<<stepNorm<<std::endl;
							std::cout<<"Residual Error: "<<resNorm<<std::endl<<std::endl;
						}
				}
			MPI_Bcast(&stepNorm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&resNorm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&option.forcing,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			sprintf(lis_option.tolerance,"-tol %f", option.forcing);
			
		}
		while(it<=option.maxIter && resNorm>option.minRes && stepNorm>option.tol);
		
		if(rank==0 && stream.verbosity==1)
			{
				fout<<"Solution: "<<std::endl;
				for (int i = 0; i < ustart.size (); ++i)
					fout << ustart[i]<<std::endl;
				
				//bim3a_print(ustart,fout);
			}
		if(rank==0 && stream.verbosity>=1)
			fout.close();

		bool converged(it<=option.maxIter);
		sol=ustart;
		return {it,resNorm,converged};
};
#endif

