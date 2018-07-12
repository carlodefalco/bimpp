#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
#include <cassert>
#include <limits>

static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return 1; }

static int
bottom_coarsening (tmesh::quadrant_iterator q)
{ return 2*(q->centroid (1) <= 0.5); }

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 5; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  tmsh.vtk_export ("p4est_interpolation_test");
  
  double delta1 = 1.5;
  double delta2 = 0.5;
  
  // Assemble matrix.
  sparse_matrix A, M;
  A.resize(tmsh.num_global_nodes());
  M.resize(tmsh.num_global_nodes());
      
  double epsilon = 1e-4;
  std::vector<double> alpha(tmsh.num_local_quadrants (), epsilon);
  std::vector<double> psi(tmsh.num_global_nodes (), 0);
      
  std::vector<double> delta(tmsh.num_local_quadrants (), 1);
  std::vector<double> zeta(tmsh.num_global_nodes (), 1);
      
  bim2a_advection_diffusion (tmsh, alpha, psi, A);
  bim2a_reaction (tmsh, delta, zeta, M);
  A += M;
      
  // Assemble right-hand side.
  std::vector<double> rhs(tmsh.num_global_nodes (), 0);
      
  std::vector<double> f(tmsh.num_local_quadrants (), 1);
  std::vector<double> g(tmsh.num_global_nodes (), 0);
      
  double x = 0, y = 0;
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii))
            {
              x = quadrant->p(0, ii);
              y = quadrant->p(1, ii);
                  
              g[quadrant->gt(ii)] = 1 - std::sinh(x / std::sqrt(epsilon)) *
                std::sinh(y / std::sqrt(epsilon)) /
                std::sinh(1 / std::sqrt(epsilon)) /
                std::sinh(1 / std::sqrt(epsilon));
            }
        }
    }
      
  // Reduce coefficients.
  std::vector<double> global_g(tmsh.num_global_nodes(), 0);
  MPI_Allreduce(g.data(), global_g.data(), g.size(),
                MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      
  bim2a_rhs (tmsh, f, global_g, rhs);
      
  // Set boundary conditions.
  func u_ex =
    [epsilon] (double x, double y)
    { return (1 - std::sinh(x / std::sqrt(epsilon)) / std::sinh(1 / std::sqrt(epsilon))) *
      (1 - std::sinh(y / std::sqrt(epsilon)) / std::sinh(1 / std::sqrt(epsilon))); };
                 
  dirichlet_bcs bcs;
  for (int i = 0; i < 4; ++i)
    bcs.push_back (std::make_tuple(0, i, u_ex));
      
  bim2a_dirichlet_bc (tmsh, bcs, A, rhs);
      
  // Solve problem.
  std::cout << "Solving linear system.";
      
  mumps mumps_solver;
      
  std::vector<double> vals;
  std::vector<int> irow, jcol;
      
  A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
      
  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);
      
  // Reduce rhs (so that rank 0 has the actual rhs).
  std::vector<double> global_rhs(tmsh.num_global_nodes(), 0);
  MPI_Reduce(rhs.data(), global_rhs.data(), rhs.size(),
             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      
  if (rank == 0)
    mumps_solver.set_rhs (global_rhs);
      
  // Solve.
  mumps_solver.analyze ();
  mumps_solver.factorize ();
  mumps_solver.solve ();
  mumps_solver.cleanup ();
      
  // Export solution.
  MPI_Bcast(global_rhs.data(), global_rhs.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  tmsh.octbin_export ("p4est_interpolation_test_u", global_rhs);
  
  std::cout << " Done." << std::endl;
  
  // Coarsen and refine.
  tmsh.set_coarsen_marker (bottom_coarsening);
  
  /*
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      tmesh::data_t * data =
	static_cast<tmesh::data_t *> (quadrant->the_quadrant->p.user_data);
         
      for (int i = 0; i < 4; ++i)
	{
	  if (i == 0)
	    {
	      std::cout << quadrant->get_global_quad_idx () << std::endl;
              for (int j = 0; j < 4; ++j)
		{
		  std::cout << data->interp_idx[j] << ", ";
		}
              std::cout << std::endl;
	    }
           
	  for (int j = 0; j < 4; ++j)
	    {
	      std::cout << data->interp_coeff[i][j] << ", ";
	    }
	  std::cout << std::endl;
	}
      std::cout << std::endl;
    }
  */
  
  tmsh.coarsen (1, partforcoarsen);
  
  // Interpolate solution at new mesh.
  std::vector<double> new_sol (tmsh.num_global_nodes ());
  interpolate_vector (tmsh, global_rhs, new_sol);
  
  MPI_Allreduce (MPI_IN_PLACE, new_sol.data (),
                 new_sol.size (), MPI_DOUBLE,
                 MPI_SUM, MPI_COMM_WORLD);
  
  tmsh.octbin_export ("p4est_interpolation_test_u_new", new_sol);
  
  MPI_Finalize ();
  
  return 0;
}
