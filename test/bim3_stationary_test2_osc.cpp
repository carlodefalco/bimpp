
/*
  Problem:	-Dnabla(u)=g
  u=1-2*x^2-2*y^2-z^2 on border
  g=6
  D=diag(0.5,0.5,1)

  Exact Solution:	u=1-2*x^2-2*y^2-z^2
*/

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>
#include <assert.h>
 
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
      std::cout << "\n\n*****\nStationary Test 2 (OSC)\n*****\n";
      
      std::cout << "read mesh" << std::endl;
      mesh msh (std::string("mesh_in_cube.msh"));

      std::cout << "export mesh" << std::endl;
      msh.write (std::string("mesh_out_cube.m"));
      
      std::cout << "compute mesh props" << std::endl;
      msh.precompute_properties ();
      
      std::cout << "assemble stiffness matrix. nnodes = " << msh.nnodes << std::endl;      
      
      bim3a_structure (msh, lhs);
      std::vector<double> dcoeff (msh.nelements*3, 1.0); //anisotropic diffusion coefficient
      for(int k=0;k<msh.nelements;++k)
        {			
          dcoeff[0+3*k]=0.5;
          dcoeff[1+3*k]=0.5;
          dcoeff[2+3*k]=1.0;
        }
      std::vector<double> ecoeff (msh.nelements, 1.0);
      std::vector<double> v (msh.nnodes, 0.0);
      std::vector<double> ncoeff (msh.nnodes, 6.0);
      
      bim3a_osc_advection_diffusion_anisotropic (msh, dcoeff, v, lhs);
      
      //bim3a_reaction (msh, ecoeff, ncoeff, lhs);
      bim3a_rhs (msh, ecoeff, ncoeff, rhs);
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
          vnodes[i]=1-2*msh.p(0,bnodes[i])*msh.p(0,bnodes[i])-2*msh.p(1,bnodes[i])*msh.p(1,bnodes[i])
            -msh.p(2,bnodes[i])*msh.p(2,bnodes[i]);
        }
      bim3a_dirichletBC(lhs,rhs,bnodes,vnodes);
			
      lhs.aij(xa,ir,jc,1);
			
      exactsolution.resize (msh.nnodes);
      for(int i=0;i<exactsolution.size();++i)
        {
          exactsolution[i]=1-2*msh.p(0,i)*msh.p(0,i)-2*msh.p(1,i)*msh.p(1,i)-msh.p(2,i)*msh.p(2,i);
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
      std::cout << "\nResult of Stationary Test \nwill be written in solution2_OSC.txt\n";
      std::ofstream fout ("solution2_OSC.txt");
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
      assert(norm < 10e-2);
    }

  mumps_solver.cleanup ();

  MPI_Finalize ();
  return (0);
}
