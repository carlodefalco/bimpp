
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
#include <fstream>
#include <stdlib.h>
//#include <mpi.h> 
LIS_INT main(LIS_INT argc, char* argv[])
{

	int int_rank,rank;
	int int_size,size;
	int n,gn,nnz,is,ie,is_r,ie_r,counter,counter_R;
  LIS_MATRIX  A;
  LIS_VECTOR  b, x;
  LIS_SOLVER  solver;
	
	LIS_INT *row,*col;
	LIS_SCALAR *value;

  LIS_INT     iter;
  double      time;
 
  lis_initialize(&argc, &argv);
	#ifdef USE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD,&int_rank);
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		rank=int_rank;
	#else
		rank=0;
		size=1;
	#endif
	
	sparse_matrix       lhs,lhs_new;
	std::vector<double> rhs1,rhs2,rhs_new;
	std::vector<int>    ir, jc;
	std::vector<double> xa;
	std::vector<double> exactsolution_start, exactsolution,sol;

	std::vector<int> bnodes;
	std::vector<double> vnodes_start, vnodes;

	std::vector<double> uold;
	double dt=1;
	int T=5;
	std::ofstream fout_sol ("Lis_Solution_TimeTest.txt");
  for (int t=1;t<=T;++t)
		{
			if (rank == 0)
		    {  
					if (t == 1)
		      	{
							std::cout << "\n\n*****\nLis test: Time Dependent Problem\n*****\n";
      
				      std::cout << "read mesh" << std::endl;
				      mesh msh (std::string("mesh_in_cube.msh"));

				      std::cout << "export mesh" << std::endl;
				      msh.write (std::string("mesh_out_cube.m"));
      
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

							std::cout << "export stiffness matrix. nnodes = " << msh.nnodes << std::endl;      
      				std::ofstream fout ("SG.m");
      				fout << lhs;
      				fout.close ();

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
			
						if(t == 1)
							lhs_new.csr(xa,jc,ir,0);
						else
							lhs_new.csr_update(xa,jc,ir,0);
						
			
						for(int i=0;i<exactsolution.size();++i)
							{
								exactsolution[i]=exactsolution_start[i];						
							}
						nnz=xa.size();
						gn=exactsolution.size();
					}
	
			#ifdef USE_MPI
				MPI_Bcast(&nnz,1,MPI_INT,0,MPI_COMM_WORLD);
				MPI_Bcast(&gn,1,MPI_INT,0,MPI_COMM_WORLD);
				
				if(rank==0)
					{
						for(int k=1;k<size;++k)
							{
								for(int i=0;i<nnz;++i)
									{
										MPI_Send(&xa[i],1,MPI_DOUBLE,k,0,MPI_COMM_WORLD);
										MPI_Send(&jc[i],1,MPI_INT,k,0,MPI_COMM_WORLD);
									}
								for(int i=0;i<gn;++i)
									{						
										MPI_Send(&ir[i],1,MPI_INT,k,0,MPI_COMM_WORLD);							
										MPI_Send(&rhs_new[i],1,MPI_DOUBLE,k,0,MPI_COMM_WORLD);
									}
								MPI_Send(&ir[gn],1,MPI_INT,k,0,MPI_COMM_WORLD);
							}
					}
				else
					{
						xa.resize(nnz);
						jc.resize(nnz);
						ir.resize(gn+1);
						rhs_new.resize(gn);
						MPI_Status *status;
		
						for(int i=0;i<nnz;++i)
							{
								MPI_Recv(&xa[i],1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,status);
								MPI_Recv(&jc[i],1,MPI_INT,0,0,MPI_COMM_WORLD,status);
							}
						for(int i=0;i<gn;++i)
							{						
								MPI_Recv(&ir[i],1,MPI_INT,0,0,MPI_COMM_WORLD,status);							
								MPI_Recv(&rhs_new[i],1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,status);
							}
						MPI_Recv(&ir[gn],1,MPI_INT,0,0,MPI_COMM_WORLD,status);
					}
			#endif
		
			if(rank==0)
				{
					n=gn/size+gn%size;
					nnz=ir[n];
					is=0;
					ie=is+nnz;
					is_r=0;
					ie_r=n+1;
				}
			else
				{
					n=gn/size;
					nnz=ir[n*(rank+1)+gn%size]-ir[n*rank+gn%size];
					is=ir[n*rank+gn%size];
					ie=is+nnz;
					is_r=gn%size+n*rank;
					ie_r=is_r+n+1;
				}
			
			row=(LIS_INT *)malloc((n+1)*sizeof(LIS_INT));
			col=(LIS_INT *)malloc(nnz*sizeof(LIS_INT));
			value=(LIS_SCALAR *)malloc(nnz*sizeof(LIS_SCALAR));
			
			//std::cout<<"rank: "<<rank<<" nnz: "<<nnz<<" n: "<<n<<" is: "<<is<<" ie: "<<ie<<" is_r: "<<is_r<<" ie_r: "<<ie_r<<" gn: "<<gn<<std::endl;
			 
			counter=0;
			counter_R=0;
			
			for(int i=is;i<ie;i++)
				{
					//std::cout<<counter<<std::endl;
					col[counter]=jc[i];
					//std::cout<<"2"<<std::endl;
					value[counter]=xa[i];
					//std::cout<<"3"<<std::endl;
					counter++;
				}
			
			for(int i=is_r;i<ie_r;++i)
				{
					if(rank!=0)
						row[counter_R]=ir[i]-ir[is_r];
					else
						row[counter_R]=ir[i];
					counter_R++;
				}
			
			lis_matrix_create(LIS_COMM_WORLD, &A);
			lis_matrix_set_size(A,n,0);
			lis_matrix_set_csr(nnz,row,col,value,A);

			//if(rank==0) std::cout<<"Assemble Matrix"<<std::endl;
			LIS_INT error=lis_matrix_assemble(A);
			if(rank==0) std::cout<<"Matrix Assembled"<<std::endl;
	
			lis_vector_create(LIS_COMM_WORLD, &b);
	
			lis_vector_set_size(b,n,0);
			lis_vector_get_range(b,&is,&ie);
			for (int i=is;i<ie;++i)
				{
					lis_vector_set_value(LIS_INS_VALUE,i,rhs_new[i],b);
				}

			lis_vector_create(LIS_COMM_WORLD, &x);
		 	
			lis_vector_duplicate(b, &x);
		 	//if(rank==0) std::cout<<"Begin solving"<<std::endl;
			lis_solver_create(&solver);
			//if(rank==0) std::cout<<"1"<<std::endl;
			lis_solver_set_optionC(solver);
			//if(rank==0) std::cout<<"2"<<std::endl;
			lis_solve(A, b, x, solver);
		 	if(rank==0) std::cout<<"System iteration solved"<<std::endl;
		
			lis_solver_get_iter(solver, &iter);
			lis_solver_get_time(solver, &time);
			 
			if(rank==0)
				{
					sol.resize(gn);
					for(int i=is;i<ie;i++)
						{
							double temp=0;
							lis_vector_get_value(x,i,&temp);
							sol[i]=temp;
						}
					MPI_Status *status;
					for(int k=1;k<size;k++)
						{				
							int is_loc,ie_loc;
							MPI_Recv(&is_loc,1,MPI_INT,k,0,MPI_COMM_WORLD,status);
							MPI_Recv(&ie_loc,1,MPI_INT,k,0,MPI_COMM_WORLD,status);
							double temp=0;
							for(int i=is_loc;i<ie_loc;i++)
								{
									MPI_Recv(&temp,1,MPI_DOUBLE,k,0,MPI_COMM_WORLD,status);
									sol[i]=temp;
								}
						}
				}
			else
				{
					MPI_Send(&is,1,MPI_INT,0,0,MPI_COMM_WORLD);
					MPI_Send(&ie,1,MPI_INT,0,0,MPI_COMM_WORLD);
					for(int i=is;i<ie;i++)
						{
							double temp=0;
							lis_vector_get_value(x,i,&temp);
							MPI_Send(&temp,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
						}

				}
		 	
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
					assert(norm < 10^-2);
				}
			//lis_vector_print(x);
			lis_solver_destroy(solver);
			lis_matrix_destroy(A);
  		lis_vector_destroy(b);
  		lis_vector_destroy(x);
		}
	fout_sol.close ();
  
  
  
	
  lis_finalize();
 
  return 0;
}

