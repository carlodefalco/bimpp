#include <bim_timing.h>
#include <tmesh_3d.h>

#include <array>
#include <cassert>
#include <cstdio>

#include "simple_connectivity_3d.h"

#include "quad_operators_3d.h"

static const int nref =  1;
static char filename[255] = "\0";

static int
uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
{ return 1; }

// my_refinement:
//
// returns 1 if the global index of node 0 is equal to 0 (i.e. it refines only
// the first quadrant) 
static int
my_refinement (tmesh_3d::quadrant_iterator quadrant)
{
  bool is_zero = false;
  if (quadrant->gt(0) == 0) 
    is_zero = true;
  return is_zero;
}

// function u
static double
my_u (double x, double y, double z)
{
  return (x*x);
}

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
  
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { tic (); }
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  if (rank == 0) { toc ("read connectivity"); }
  MPI_Barrier (MPI_COMM_WORLD);

  // definition of u
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

  // print u
  if (rank == 0)
    {
      std::cout << "\n\n***** u: " << std::endl;
      for (double d : u_vec)
        std::cout << d << " ";
      std::cout << std::endl;
    }

  // computation of the gradient
  gradient3 grad_vec = bim2c_quadtree_pde_recovered_gradient(tmsh,u_vec);
  q1_vec dudx_1 = std::get<0>(grad_vec);
  q1_vec dudy_1 = std::get<1>(grad_vec);  
  q1_vec dudz_1 = std::get<2>(grad_vec);

  // print the gradient
  if (rank == 0)
  {
    std::cout << "\n\n***** dudx_1: " << std::endl;
    for (double d : dudx_1)
      std::cout << d << " ";
    std::cout << std::endl;
    std::cout << "\n\n***** dudy_1: " << std::endl;
    for (double d : dudy_1)
      std::cout << d << " ";
    std::cout << std::endl;
    std::cout << "\n\n***** dudz_1: " << std::endl;
    for (double d : dudz_1)
      std::cout << d << " ";
    std::cout << std::endl;
  }

  // estimator gradient
  std::vector<double> est_grad;
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
      est_grad.push_back(estimator_grad(q, grad_vec, u_vec));

  // print estimator gradient  
  if (rank == 0)
    {
      std::cout << "\n***** estimator grad:" << std::endl;
      for (double d : est_grad)
        std::cout << d << " ";
      std::cout << std::endl;
    }  

  // first mesh uniform refinement
  for (auto i = 0; i < 1; ++i)
    {

      MPI_Barrier (MPI_COMM_WORLD);
      if (rank == 0) { tic (); }

      tmsh.set_refine_marker (uniform_refinement);
      recursive = 0; partforcoarsen = 1;
      tmsh.refine (recursive, partforcoarsen);

      if (rank == 0) { toc ("uniform refinement"); }
      MPI_Barrier (MPI_COMM_WORLD);

      sprintf (filename, "myexe_initial_mesh_%5.5d", i);
      tmsh.vtk_export (filename);
    }
/*
  // mesh info
  if (rank == 0)
  	{
  	  std::cout << "**************** mesh info ****************" << std::endl;
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
	          std:: cout << "\tnode: " << nn 
	          			 << " x: " << x
	          			 << " y: " << y
	          			 << " z: " << z
	          			 << " gt: " << gidx << std::endl;
	        }  
	      std::cout << std::endl; 
	      ++nq; 
	    } 
	} 	  
*/
  // update of u
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

  // print u
  if (rank == 0)
    {
      std::cout << "\n\n***** u: " << std::endl;
      for (double d : u_vec)
        std::cout << d << " ";
      std::cout << std::endl;
    }    

  // computation of the gradient with the intermediate mesh
  gradient3 grad_vec_2 = bim2c_quadtree_pde_recovered_gradient(tmsh,u_vec);
  q1_vec dudx_2 = std::get<0>(grad_vec_2);
  q1_vec dudy_2 = std::get<1>(grad_vec_2);  
  q1_vec dudz_2 = std::get<2>(grad_vec_2);
  // print the gradient
  if (rank == 0)
  {
    std::cout << "\n\n***** dudx_2: " << std::endl;
    for (double d : dudx_2)
      std::cout << d << " ";
    std::cout << std::endl;
    std::cout << "\n\n***** dudy_2: " << std::endl;
    for (double d : dudy_2)
      std::cout << d << " ";
    std::cout << std::endl;
    std::cout << "\n\n***** dudz_2: " << std::endl;
    for (double d : dudz_2)
      std::cout << d << " ";
    std::cout << std::endl;
  }

  // estimator gradient
  std::vector<double> est_grad_2;
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
  	est_grad_2.push_back(estimator_grad(q, grad_vec_2, u_vec));

  // print estimator gradient  
  if (rank == 0)
    {
      std::cout << "\n***** estimator grad_2:" << std::endl;
      for (double d : est_grad_2)
        std::cout << d << " ";
      std::cout << std::endl;
    }     

  // second mesh refinement
  for (auto i = 0; i < nref; ++i)
    {
      MPI_Barrier (MPI_COMM_WORLD);
      if (rank == 0) { tic (); }

      tmsh.set_refine_marker (my_refinement);
      recursive = 1; partforcoarsen = 1;
      tmsh.refine (recursive, partforcoarsen);

      if (rank == 0) { toc ("my refinement"); }
      MPI_Barrier (MPI_COMM_WORLD);

      sprintf (filename, "myexe_%5.5d", i);
      tmsh.vtk_export (filename);
    } 

  // update of u
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

  // print u
  if (rank == 0)
    {
      std::cout << "\n\n***** u: " << std::endl;
      for (double d : u_vec)
        std::cout << d << " ";
      std::cout << std::endl;
    }      

  // computation of the gradient with the final mesh
  gradient3 grad_vec_3 = bim2c_quadtree_pde_recovered_gradient(tmsh,u_vec);
  q1_vec dudx_3 = std::get<0>(grad_vec_3);
  q1_vec dudy_3 = std::get<1>(grad_vec_3);  
  q1_vec dudz_3 = std::get<2>(grad_vec_3);
  // print the gradient
  if (rank == 0)
  {
    std::cout << "\n\n***** dudx_3: " << std::endl;
    for (double d : dudx_3)
      std::cout << d << " ";
    std::cout << std::endl;
    std::cout << "\n\n***** dudy_3: " << std::endl;
    for (double d : dudy_3)
      std::cout << d << " ";
    std::cout << std::endl;
    std::cout << "\n\n***** dudz_3: " << std::endl;
    for (double d : dudz_3)
      std::cout << d << " ";
    std::cout << std::endl;
  }

  // estimator gradient
  std::vector<double> est_grad_3;
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
      est_grad_3.push_back(estimator_grad(q, grad_vec_3, u_vec));

  // print estimator gradient  
  if (rank == 0)
    {
      std::cout << "***** estimator grad_3:" << std::endl;
      for (double d : est_grad_3)
        std::cout << d << " ";
      std::cout << std::endl;
    }  

  if (rank == 0) {print_timing_report();}
  MPI_Barrier (MPI_COMM_WORLD);
  
  MPI_Finalize ();
  return 0;

}

