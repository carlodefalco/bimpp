
/*
  Problem:	-nabla(u)+u=g
  u=sin(pi*x)+sin(pi*y)+sin(pi*z) on border
  g=(pi^2+1)*(sin(pi*x)+sin(pi*y)+sin(pi*z))

  Exact Solution:	u=sin(pi*x)+sin(pi*y)+sin(pi*z)
*/

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <bim_config.h>

int main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);    

  sparse_matrix       lhs;
  std::vector<double> rhs;
  std::vector<int>    ir, jc;
  std::vector<double> xa;
  std::vector<double> exactsolution;
  if (rank == 0)
    {  
      std::cout << "\n\n*****\nStationary Test 4\n*****\n";
			      
      std::cout << "read mesh" << std::endl;
      mesh msh (data_dir + std::string("mesh_in_cube2.msh"));

      std::cout << "export mesh" << std::endl;
      msh.write (std::string("mesh_out_cube2.m"));
      
      std::cout << "compute mesh props" << std::endl;
      msh.precompute_properties ();
      
      std::cout << "assemble stiffness matrix. nnodes = " << msh.nnodes << std::endl;      
      
      bim3a_structure (msh, lhs);
      std::vector<double> ecoeff (msh.nelements, 1.0); //isotropic diffusion coefficient
      std::vector<double> v (msh.nnodes, 0.0);
      std::vector<double> ncoeff (msh.nnodes, 1.0);
      std::vector<double> nodecoeff (msh.nnodes, 0.0);      

      for (int i=0; i< msh.nnodes;++i)
        nodecoeff[i]=(M_PI*M_PI+1)*(sin(M_PI*msh.p(0,i))+sin(M_PI*msh.p(1,i))+sin(M_PI*msh.p(2,i)));

      bim3a_advection_diffusion (msh, ecoeff, v, lhs);
      
      bim3a_reaction (msh, ecoeff, ncoeff, lhs);
      bim3a_rhs (msh, ecoeff, nodecoeff, rhs);
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
			
      std::vector<int> bnodes;
			
      std::vector<double> vnodes;
      bim3a_boundary_nodes(msh,sidelist,bnodes);
			
      vnodes.resize(bnodes.size());
      for(int i=0;i<vnodes.size();++i)
        {
          vnodes[i]=sin(M_PI*msh.p(0,bnodes[i]))+sin(M_PI*msh.p(1,bnodes[i]))+sin(M_PI*msh.p(2,bnodes[i]));
								
        }
      bim3a_dirichletBC(lhs,rhs,bnodes,vnodes);
			
      lhs.aij(xa,ir,jc,1);
			
      exactsolution.resize (msh.nnodes);
      for(int i=0;i<exactsolution.size();++i)
        {
          exactsolution[i]=sin(M_PI*msh.p(0,i))+sin(M_PI*msh.p(1,i))+sin(M_PI*msh.p(2,i));
        }
    }

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
    {
      std::cout << "\nResult of Stationary Test \nwill be written in solution4.txt\n";
      std::ofstream fout ("solution4.txt");
      fout << std::endl;
			
      double norm=0;
      //double normexact=0;
      for (int k = 0; k < rhs.size (); ++k)
        {
          fout << rhs[k] << "  " << exactsolution[k]<< std::endl;
          norm+=(exactsolution[k]-rhs[k])*(exactsolution[k]-rhs[k]);
          //normexact+=(exactsolution[k])*(exactsolution[k]);
        }
      fout.close ();
      std::cout<<"Error: "<<norm<<std::endl;
      assert(norm < 10e-1);
    }

  mumps_solver.cleanup ();

  MPI_Finalize ();
  return (0);
}
