/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/

/*
	Problem:	du/dt-div(D(grad(u)-grad(v)*u))=g
						u=(1-x^2-y^2-z^2)*exp(-t) on border
						g=(11+x^2+y^2+z^2-2x-2y-2z)*exp(-t)
						D=diag(1.0,2.0,3.0)
						v=[1.0,0.5,1/3]

	Exact Solution:	u=(1-x^2-y^2-z^2)*exp(-t)
*/

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>
#include <cmath>
  
int main (void)
{

  MPI::Init ();
  int rank = MPI::COMM_WORLD.Get_rank ();
  int size = MPI::COMM_WORLD.Get_size ();
    
	sparse_matrix       lhs;
	sparse_matrix				lhs_new;
	std::vector<double> rhs1;
	std::vector<double> rhs2;
	std::vector<double> rhs_new;
	std::vector<int>    ir, jc;
	std::vector<double> xa;

  std::vector<double> exactsolution_start, exactsolution;

	std::vector<int> bnodes;
	std::vector<double> vnodes_start, vnodes;

	std::vector<double> uold;
	double dt=0.2;
	int T=5;
  for (int t=1;t<=T;++t)
		{
			if (rank == 0)
		    {  
					if (t == 1)
		      	{
							std::cout << "\n\n*****\nTime Dependent Test 3\n*****\n";
      
				      std::cout << "read mesh" << std::endl;
				      mesh msh (std::string("mesh_in_cube.msh"));

				      std::cout << "export mesh" << std::endl;
				      msh.write (std::string("mesh_out_cube.m"));
      
				      std::cout << "compute mesh props" << std::endl;
				      msh.precompute_properties ();
      
				      std::cout << "assemble stiffness matrix. nnodes = " << msh.nnodes << std::endl;      
      
							bim3a_structure (msh, lhs);
				      std::vector<double> ecoeff (msh.nelements, 1.0); 
							std::vector<double> dcoeff (msh.nelements*3,1.0);//anisotropic diffusion coefficients
							for (int i=0; i < msh.nelements;++i)
								{
									dcoeff[0+3*i]=1.0;
									dcoeff[1+3*i]=2.0;
									dcoeff[2+3*i]=3.0;
								}
				      std::vector<double> v (msh.nnodes, 0.0);
							for (int i=0; i < msh.nnodes;++i)
								v[i]=msh.p(0,i)+0.5*msh.p(1,i)+msh.p(2,i)/3.0;

							std::vector<double> ncoeff(msh.nnodes, 1/dt);
							std::vector<double> nodecoeff1(msh.nnodes,1/dt);
							std::vector<double> nodecoeff2(msh.nnodes,0.0);

							for(int i=0;i<msh.nnodes;++i)
								nodecoeff2[i]=12-1+msh.p(0,i)*msh.p(0,i)+msh.p(1,i)*msh.p(1,i)+msh.p(2,i)*msh.p(2,i)
															-2*msh.p(0,i)-2*msh.p(1,i)-2.0*msh.p(2,i);

							bim3a_advection_diffusion_anisotropic (msh, dcoeff, v, lhs);
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
						}
								
						rhs_new.resize(rhs1.size());
						lhs_new=lhs;

						for (int i=0;i<rhs_new.size();++i)
							{
								rhs_new[i]=rhs2[i]*exp(-t*dt)+rhs1[i]*uold[i];
							}
						for (int i=0;i<vnodes.size();++i)
							{	
								vnodes[i]=vnodes_start[i]*exp(-t*dt);
							}
								
						bim3a_dirichletBC(lhs_new,rhs_new,bnodes,vnodes);
			
						lhs_new.aij(xa,ir,jc,1);
			
						for (int i=0;i<exactsolution.size();++i)
							{
								exactsolution[i]=exactsolution_start[i]*exp(-t*dt);							
							}
					}	
    		

  		mumps mumps_solver;
  
  		if (rank == 0)
				mumps_solver.set_lhs_structure (lhs_new.rows (), ir, jc);
  
  		mumps_solver.analyze ();
  
  		if (rank == 0)
    		mumps_solver.set_lhs_data (xa);

  		mumps_solver.factorize ();
  
  		if (rank == 0)
				mumps_solver.set_rhs (rhs_new);

  		mumps_solver.solve ();
  
  		if (rank == 0)
    		{
					char nomefile[50];
					sprintf(nomefile,"Solution_TD3_%02d",t);
					if(t==1) 
						std::cout <<"\nResult of Time Dependent Test\n\nIteration 1 will be written in " << nomefile <<".txt\n";
					else 
						std::cout << "\nIteration " << t <<" will be written in " << nomefile <<".txt\n";
      		std::ofstream fout (nomefile);
      		fout << std::endl;
					double norm=0;
					//double normexact=0;
      		for (int k = 0; k < rhs_new.size (); ++k)
						{
	      			fout << rhs_new[k] << "  " << exactsolution[k]<< std::endl;
							norm+=(exactsolution[k]-rhs_new[k])*(exactsolution[k]-rhs_new[k]);
							//normexact+=(exactsolution[k])*(exactsolution[k]);
							uold[k]=rhs_new[k];
						}
      		fout.close ();
					std::cout<<"Error: "<<norm<<std::endl;
    		}

			mumps_solver.cleanup ();

		}
	MPI::Finalize();
  return (0);
}
