
/*
  Problem:      du/dt-nabla(u)=g
  u=1-x^2-y^2-z^2 on border
  g=6

  Exact Solution:       u=1-x^2-y^2-z^2
*/

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <bim_config.h>

int main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  
  sparse_matrix       lhs;
  sparse_matrix       lhs_new;
  std::vector<double> rhs1;
  std::vector<double> rhs2;
  std::vector<double> rhs_new;
  std::vector<int>    ir, jc;
  std::vector<double> xa;
  std::vector<double> exactsolution_start, exactsolution;

  std::vector<int> bnodes;
  std::vector<double> vnodes_start, vnodes;

  std::vector<double> uold;
  double dt=1;
  int T=5;
  char nomefile[]="Solution_TimeTest1.txt";
  std::ofstream fout (nomefile);
  for (int t=1;t<=T;++t)
    {
      if (rank == 0)
        {  
          if (t == 1)
            {
              std::cout << "\n\n*****\nTime Dependent Test\n*****\n";
      
              std::cout << "read mesh" << std::endl;
              mesh msh (data_dir + std::string ("mesh_in_cube.msh"));

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
            }
                                                                
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
            lhs_new.aij(xa,ir,jc,1);
          else
            lhs_new.aij_update(xa,ir,jc,1);
                        
          for(int i=0;i<exactsolution.size();++i)
            {
              exactsolution[i]=exactsolution_start[i];                                          
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
          //sprintf(nomefile,"Solution_TD1_%02d.txt",t);
          if(t==1) 
            std::cout <<"\nResult of Time Dependent Test\nwill be written in " << nomefile <<std::endl;
                                         
          std::cout << "\nIteration " << t <<std::endl;
                
          fout << "Iteration " << t <<std::endl;
          double norm=0;
          //double normexact=0;
          for (int k = 0; k < rhs_new.size (); ++k)
            {
              fout << rhs_new[k] << "  " << exactsolution[k]<< std::endl;
              norm+=(exactsolution[k]-rhs_new[k])*(exactsolution[k]-rhs_new[k]);
              //normexact+=(exactsolution[k])*(exactsolution[k]);
              uold[k]=rhs_new[k];
            }
          fout << std::endl;
          std::cout<<"Error: "<<norm<<std::endl;
        }
      mumps_solver.cleanup ();
    }
  fout.close ();
  MPI_Finalize ();
  return (0);
}
