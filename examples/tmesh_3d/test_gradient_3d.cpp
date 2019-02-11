#include <bim_timing.h>
#include <tmesh_3d.h>

#include <array>
#include <cassert>
#include <cstdio>
#include <string>

#include "simple_connectivity_3d.h"   // one tree

#include "quad_operators_3d.h"

static const int nref =  1;
static char filename[255] = "\0";

// uniform_refinement:
//
// returns 1 ----> all quadrants are refined
static int
uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
{ return 1; }

// my_refinement:
//
// returns 1 if the global index of node 0 is equal to 0 ----> only the first 
// quadrant is refined
static int
my_refinement (tmesh_3d::quadrant_iterator quadrant)
{
  bool is_zero = false;
  if (quadrant->gt(0) == 0) 
    is_zero = true;
  return is_zero;
}

// my_u:
//
// implements the function u
static double
my_u (double x, double y, double z)
{
  return (5*x + 4*y + 2*z);
}

// print_mesh_info:
//
// it prints the number of nodes and quadrants of the mesh; it scans all the 
// quadrants and prints all the information about the nodes (: local index,
// global index, coordinates)
static void
print_mesh_info (tmesh_3d& tmsh, const std::string& str)
{
  std::cout << "**************** mesh info " << str << " ****************" 
            << std::endl;
  std::cout << "num owned nodes " << tmsh.num_owned_nodes() << std::endl;
  std::cout << "num local nodes " << tmsh.num_local_nodes() << std::endl;
  std::cout << "num global nodes " << tmsh.num_global_nodes() << std::endl;
  std::cout << "num local quads " << tmsh.num_local_quadrants() 
            << std::endl;
  std::cout << "num global quads " << tmsh.num_global_quadrants() 
            << std::endl;
	  unsigned nq = 0;
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
    {
      std::cout << "quadrant " << nq << std::endl;	
      for (int nn = 0; nn < 8; ++nn)
        {
          double x = q->p(0,nn);
          double y = q->p(1,nn);
          double z = q->p(2,nn);
          unsigned gidx = q->gt(nn);
          std::cout << "\tnode: " << nn 
          			 << " x: " << x
          			 << " y: " << y
          			 << " z: " << z
          			 << " gt: " << gidx << std::endl;
        }  
      std::cout << std::endl; 
      ++nq; 
    }
}

// print:
//
// prints on std::cout an object of type TT (a vector or an array of doubles) 
template <typename TT>
void print (const TT& vec, const std::string& str)
{
  std::cout << "\n***** " << str << ": " << std::endl;
  for (double d : vec)
    std::cout << d << " ";
  std::cout << std::endl << std::endl;
}

// MPI_User_function.
static void replace(double *invec, double *inoutvec,
                    int *len, MPI_Datatype *dtype)
{
  for (int i = 0; i < *len; ++i)
    if (/*invec[i] != 0 && inoutvec[i] == 0*/ i < 15)
      inoutvec[i] = invec[i];
}



/*********************************** MAIN ***********************************/
int
main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  // set to true if you want the information about the mesh to be printed 
  // after any refinement
  bool print_mesh = false;


  /************************ initialization of the mesh ***********************/

  MPI_Barrier (mpicomm);
  if (rank == 0) { tic (); }
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  if (rank == 0) { toc ("read connectivity"); }
  MPI_Barrier (mpicomm);

  // definition of u (initial mesh)
  q1_vec u_vec(tmsh.num_global_nodes());
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
    {
      for (int nn = 0; nn < 8; ++nn)
        {
          if (! q->is_hanging(nn))
          {
            double x = q->p(0,nn);
            double y = q->p(1,nn);
            double z = q->p(2,nn);
            double uu = my_u(x,y,z);
            u_vec[q->gt(nn)] = uu;
          }
        }  
    }

  // print u (initial mesh)
  MPI_Barrier (mpicomm);
  if (rank == 0)
    print(u_vec,"u");   



  /*********************** first mesh uniform refinement *********************/
 
  for (auto i = 0; i < nref; ++i)
    {

      MPI_Barrier (mpicomm);
      if (rank == 0) { tic (); }

      tmsh.set_refine_marker (uniform_refinement);
      recursive = 0; partforcoarsen = 1;
      tmsh.refine (recursive, partforcoarsen);

      if (rank == 0) { toc ("uniform refinement"); }
      MPI_Barrier (mpicomm);

      sprintf (filename, "test_grad_3d_initial_mesh_%5.5d", i);
      tmsh.vtk_export (filename);
    }

  // mesh info (intermediate mesh - uniform refinement)
  if (print_mesh)
	{  
	  MPI_Barrier (mpicomm);
	  if (rank == 0)
	  	print_mesh_info(tmsh,"after uniform refinement");  
	}	  

  // update of u (intermediate mesh - uniform refinement)
  u_vec.resize(tmsh.num_global_nodes());
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
    {
      for (int nn = 0; nn < 8; ++nn)
        {
          if (! q->is_hanging(nn))
            {
              double x = q->p(0,nn);
              double y = q->p(1,nn);
              double z = q->p(2,nn);
              double uu = my_u(x,y,z);
              u_vec[q->gt(nn)] = uu;
            }
        }  
    }

  // print u (intermediate mesh - uniform refinement)
  MPI_Barrier (mpicomm);  
  if (rank == 0)
    print(u_vec,"u");

  // computation of the gradient (intermediate mesh - uniform refinement)
  gradient3 grad_vec_2 = bim2c_quadtree_pde_recovered_gradient(tmsh,u_vec);
  q1_vec dudx_2 = std::get<0>(grad_vec_2);
  q1_vec dudy_2 = std::get<1>(grad_vec_2);  
  q1_vec dudz_2 = std::get<2>(grad_vec_2);
  
  // print the gradient (intermediate mesh - uniform refinement)
  MPI_Barrier (mpicomm);
  if (rank == 0)
  {
    print(dudx_2,"dudx_2");
    print(dudy_2,"dudy_2");
    print(dudz_2,"dudz_2");
  }

  // estimator gradient (intermediate mesh - uniform refinement)
  q1_vec est_grad_2;
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
  	est_grad_2.push_back(estimator_grad(q, grad_vec_2, u_vec));

  // print estimator gradient (intermediate mesh - uniform refinement)
  MPI_Barrier (mpicomm);
  if (rank == 0)
    print(est_grad_2,"estimator grad _2");

  // recovered solution (intermediate mesh - uniform refinement)
  q2_vec3 u_rec_2 = bim2c_quadtree_pde_recovered_solution (tmsh,u_vec,
  														   grad_vec_2);  

  // print recovered solution (intermediate mesh - uniform refinement)
  MPI_Barrier (mpicomm);
  if (rank == 0)
  	{
  	  std::cout << "\n***** recovered solution _2:" << std::endl;
  	  for (const auto& arr : u_rec_2)
  	  	print(arr,"");	
  	}

  // estimator solution (intermediate mesh - uniform refinement)
  q1_vec est_sol_2;
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
  	est_sol_2.push_back(estimator_sol(q, u_rec_2, u_vec));

  // print estimator solution (intermediate mesh - uniform refinement)  
  MPI_Barrier (mpicomm);
  if (rank == 0)
    print(est_sol_2,"estimator sol _2");




  /************************** second mesh refinement *************************/

  for (auto i = 0; i < nref; ++i)
    {
      MPI_Barrier (mpicomm);
      if (rank == 0) { tic (); }

      tmsh.set_refine_marker (my_refinement);
      recursive = 1; partforcoarsen = 1;
      tmsh.refine (recursive, partforcoarsen);

      if (rank == 0) { toc ("my refinement"); }
      MPI_Barrier (mpicomm);

      sprintf (filename, "test_grad_3d_%5.5d", i);
      tmsh.vtk_export (filename);
    } 

  // mesh info (final mesh - my_refinement)
  if (print_mesh)
	 {
  	  MPI_Barrier (mpicomm);
  	  if (rank == 0)
  	    print_mesh_info(tmsh,"after my_refinement");
    } 

  // update of u (final mesh - my_refinement)
  u_vec.resize(tmsh.num_global_nodes());
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
    {
      for (int nn = 0; nn < 8; ++nn)
        {
          if (! q->is_hanging(nn))
            {
              double x = q->p(0,nn);
              double y = q->p(1,nn);
              double z = q->p(2,nn);
              double uu = my_u(x,y,z);
              u_vec[q->gt(nn)] = uu;
            }
        }  
    } 

  // Send u data to all processes so that non-assigned values
  // on current rank get assigned by other ranks.
  MPI_Op op;
  MPI_Op_create ((MPI_User_function *) replace, 1, &op);

  MPI_Allreduce (MPI_IN_PLACE, u_vec.data (), u_vec.size (), MPI_DOUBLE,
                 op, mpicomm);  

  // print u (final mesh - my_refinement)
  MPI_Barrier (mpicomm);
  if (rank == 0)
    print(u_vec,"u");      

  // computation of the gradient (final mesh - my_refinement)
  gradient3 grad_vec_3 = bim2c_quadtree_pde_recovered_gradient(tmsh,u_vec);
  q1_vec dudx_3 = std::get<0>(grad_vec_3);
  q1_vec dudy_3 = std::get<1>(grad_vec_3);  
  q1_vec dudz_3 = std::get<2>(grad_vec_3);
  
  // print the gradient (final mesh - my_refinement)
  MPI_Barrier (mpicomm);
  if (rank == 0)
  {
    print(dudx_3,"dudx_3");
    print(dudy_3,"dudy_3");
    print(dudz_3,"dudz_3");
  }

  // estimator gradient (final mesh - my_refinement)
  q1_vec est_grad_3;
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
      est_grad_3.push_back(estimator_grad(q, grad_vec_3, u_vec));

  // print estimator gradient (final mesh - my_refinement)
  MPI_Barrier (mpicomm);
  if (rank == 0)
    print(est_grad_3,"estimator grad _3");  

  // recovered solution (final mesh - my_refinement)
  q2_vec3 u_rec_3 = bim2c_quadtree_pde_recovered_solution (tmsh,u_vec,
  														   grad_vec_3);  

  // print recovered solution (final mesh - my_refinement)
  MPI_Barrier (mpicomm);
  if (rank == 0)
  	{
  	  std::cout << "\n***** recovered solution _3:" << std::endl;
  	  for (const auto& arr : u_rec_3)
  	  	print(arr,"");	
  	}

  // estimator solution (final mesh - my_refinement)
  q1_vec est_sol_3;
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
  	est_sol_3.push_back(estimator_sol(q, u_rec_3, u_vec));

  // print estimator solution (final mesh - my_refinement) 
  MPI_Barrier (mpicomm);
  if (rank == 0)
    print(est_sol_3,"estimator sol _3");

  if (rank == 0) {print_timing_report();}
  MPI_Barrier (mpicomm);
  
  MPI_Finalize ();
  return 0;

}