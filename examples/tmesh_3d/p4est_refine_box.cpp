#include <quad_operators_3d.h>

#include <bim_timing.h>

#include <cstdio>
#include <algorithm>
#include <cmath>

char name[255];
char step[255];

static int
uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
{ return 1; }

static int
rectangle_list_refinement (tmesh_3d::quadrant_iterator quadrant,
                           double L, double h)
{
  double x0 = quadrant->p (0, 0);
  double y0 = quadrant->p (1, 0);
  double z0 = quadrant->p (2, 0);

  double x1 = quadrant->p (0, 7);
  double y1 = quadrant->p (1, 7);
  double z1 = quadrant->p (2, 7);

  double l = .5 * (L - 3. * h);
  
  if (x0>h & x1<h+l
      & y0>h & y1<h+l
      & z0>h & z1<h+l)
    return true;
  else
    return false;
}


int main(int argc, char ** argv)
{
  MPI_Init (&argc, &argv);
  
  int      recursive, partforcoarsen, balance;
  MPI_Comm mpicomm = MPI_COMM_WORLD;  							
  int      rank, size;
  tmesh_3d tmsh;

  size_t num_quads = 0;			// number of quadrants
  size_t nq = 0;						// index of current quadrant
  
  double local_h = 0;				// mesh size on current processor
  double global_h = 0;			// global mesh size
  
  using q1_vec = q1_vec<distributed_vector>;
  using gradient3 = gradient3<distributed_vector>;
  using idx_t = tmesh_3d::idx_t;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);
    
  MPI_Barrier (MPI_COMM_WORLD); { if (rank == 0) tic (); }
  
  // Create mesh.
  /*
  std::vector<double> p = {0, 0, 0,
                           0, 0, 1,
                           0, 1, 0,
                           0, 1, 1,
                           1, 0, 0,
                           1, 0, 1,
                           1, 1, 0,
                           1, 1, 1};
  */
  std::vector<double> p = {0, 0, 0,
                           1, 0, 0,
                           0, 1, 0,
                           1, 1, 0,
                           0, 0, 1,
                           1, 0, 1,
                           0, 1, 1,
                           1, 1, 1};
  std::vector<int> t = {1, 2, 3, 4, 5, 6, 7, 8, 1};
  
  tmsh.read_connectivity (&(p[0]), 8, &(t[0]), 1);

  double L = 1., h = 1./5.;
  
  std::function<int (tmesh_3d::quadrant_iterator)> box_refinement =
    [L,h] (tmesh_3d::quadrant_iterator qi)
    { return rectangle_list_refinement(qi, L, h); };

  if (rank == 0) { toc ("*** Initialization ***"); }

  // Definition of exact solution
  func3 u_ex = 
    [L] (double x, double y, double z) -> double
    {
      return (sin(2*M_PI*x/L) * cos(2*M_PI*y/L) * sin(2*M_PI*z/L));
      //return (5*x + 4*y + 2*z);
    };

  // Definition of exact derivatives
  func3 du_x_ex =
    [L] (double x, double y, double z) -> double
    {
      return (2*M_PI * cos(2*M_PI*x/L) * cos(2*M_PI*y/L) * sin(2*M_PI*z/L) / L);
      //return 5.;
    };

  func3 du_y_ex =
    [L] (double x, double y, double z) -> double
    {
      return (-2*M_PI * sin(2*M_PI*x/L) * sin(2*M_PI*y/L) * sin(2*M_PI*z/L) /L);
      //return 4.;
    };

  func3 du_z_ex =
    [L] (double x, double y, double z) -> double
    {
      return (2*M_PI * sin(2*M_PI*x/L) * cos(2*M_PI*y/L) * cos(2*M_PI*z/L) / L);
      //return 2.;
    };

  // Containers for errors
  std::array<double,4> loc_err;			// loc_err[0] -> error
  																	// loc_err[1] -> error_du
  																	// loc_err[2] -> error_star
    																// loc_err[3] -> error_star_du

  std::array<double,4> errors;			// errors[0] -> errors
  																	// errors[1] -> errors_du
  																	// errors[2] -> errors_star
    																// errors[3] -> errors_star_du

  std::vector<idx_t> nnodes;
  std::vector<idx_t> nquads;
  std::vector<double> h_vec;

  // Uniform refinement.
  recursive = 0;
  partforcoarsen = 1;
  
  int file_number = 0, refine_number = 0;
  for (int cycle = 0; cycle < 1; ++cycle)
    {
        MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
        tmsh.set_refine_marker (uniform_refinement);
        tmsh.refine(recursive, partforcoarsen);        
        if (rank == 0) 
          { 
             sprintf (step, "*** Refinement and balancing %3.3d ***", 
                      refine_number++);
             toc (step); 
          }

        MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
        sprintf(name,"p4est_refine_box_%4.4d",file_number++);
        tmsh.vtk_export (name);
        if (rank == 0) { toc ("*** Export ***"); }
    }

  // *************************************************************************

  // definition of u (uniform refined mesh)
  q1_vec u_vec(tmsh.num_owned_nodes());
  bim3a_solution_with_ghosts(tmsh,u_vec);
  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
    {
      for (int nn = 0; nn < 8; ++nn)
        {
          // assemble non-hanging nodes
          if (! q->is_hanging(nn))
            u_vec[q->gt(nn)] = u_ex(q->p(0,nn), q->p(1,nn), q->p(2,nn));
        }  
    }
  u_vec.assemble(replace_op);

  // Compute gradient (uniform refined mesh)
  gradient3 grad_star = bim3c_quadtree_pde_recovered_gradient(tmsh,u_vec);

  // Compute recovered solution (uniform refined mesh)
  q2_vec3 u_star = bim3c_quadtree_pde_recovered_solution(tmsh,u_vec,
  																												grad_star);

  // Compute h and error norms (uniform refined mesh)
  num_quads = tmsh.num_global_quadrants();
  local_h = std::numeric_limits<double>::max ();
  global_h = 0;
  loc_err.fill(0.);
  errors.fill(0.);

  for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
    {
    	double hx = q->p(0, 7) - q->p(0, 0);
      double hy = q->p(1, 7) - q->p(1, 0);
      double hz = q->p(2, 7) - q->p(2, 0);
      local_h = std::min(local_h, std::sqrt(hx*hx + hy*hy + hz*hz));

      // || u_k - u_ex ||_{L^2}
      loc_err[0] += std::pow(l2_error(q, u_ex, u_vec), 2);
      // || u^*_k - u_ex ||_{L^2}
      loc_err[1] += std::pow(l2_star_error(q, u_ex, u_star), 2);
      // | sigma - grad(u_ex) |_{H^1}
      loc_err[2] += std::pow(semih1_error(q, du_x_ex, du_y_ex, du_z_ex, 
      																		u_vec), 2);
      // | sigma^* - grad(u_ex) |_{H^1}
      loc_err[3] += std::pow(semih1_star_error(q, du_x_ex,du_y_ex,du_z_ex, 
      																					grad_star), 2);
    }

  MPI_Reduce(&local_h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);

  MPI_Reduce(loc_err.data(), errors.data(), 4, MPI_DOUBLE, MPI_SUM,
  						0, mpicomm);

  for (double& d : errors)
  	d = std::sqrt(d);

  nnodes.push_back(tmsh.num_global_nodes());
  nquads.push_back(num_quads);
  h_vec.push_back(global_h);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    {
    	std::cout << "nnodes = " << nnodes[nnodes.size()-1] << std::endl;
    	std::cout << "nquads = " << nquads[nquads.size()-1] << std::endl;
    	std::cout << "h = " << h_vec[h_vec.size()-1] << std::endl;
      for (unsigned j=0; j<4; ++j)
        std::cout << "err" << j << " = " << errors[j] << std::endl;
      std::cout << std::endl;
    }

  // *************************************************************************  
  
  // Refine according to boxes.  
  for (int cycle = 0; cycle < 4; ++cycle)
    {
      // Adaptive refinement.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }

      tmsh.set_refine_marker (uniform_refinement); 
      tmsh.refine (recursive, partforcoarsen);
      
      tmsh.set_refine_marker (box_refinement);
      tmsh.refine (recursive, partforcoarsen);
      
      if (rank == 0) 
          { 
             sprintf (step, "*** Refinement and balancing %3.3d ***", 
                      refine_number++);
             toc (step); 
          }

        // *******************************************************************

        // definition of u (adaptive refined mesh)
        q1_vec u_vec_1(tmsh.num_owned_nodes());
        bim3a_solution_with_ghosts(tmsh,u_vec_1);
        for (auto q = tmsh.begin_quadrant_sweep();
                  q != tmsh.end_quadrant_sweep();
                  ++q)
          {
            for (int nn = 0; nn < 8; ++nn)
              {
                // assemble non-hanging nodes
                if (! q->is_hanging(nn))
                  u_vec_1[q->gt(nn)] = u_ex(q->p(0,nn), q->p(1,nn), 
                                            q->p(2,nn));
              }  
          }
        u_vec_1.assemble(replace_op);

        // gradient (adaptive refined mesh)
        gradient3 grad_star_1 = bim3c_quadtree_pde_recovered_gradient(tmsh,
                                                                      u_vec_1);

        // recovered solution (adaptive refined mesh)
        q2_vec3 u_star_1 = bim3c_quadtree_pde_recovered_solution(tmsh,u_vec_1,
                                                                  grad_star_1);

        // Compute h and error norms (adaptive refined mesh)
			  num_quads = tmsh.num_global_quadrants();
			  local_h = std::numeric_limits<double>::max ();
			  global_h = 0;
			  loc_err.fill(0.);
			  errors.fill(0.);

        for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
			    {
			    	double hx = q->p(0, 7) - q->p(0, 0);
			      double hy = q->p(1, 7) - q->p(1, 0);
			      double hz = q->p(2, 7) - q->p(2, 0);
			      local_h = std::min(local_h, std::sqrt(hx*hx + hy*hy + hz*hz));

			      // || u_k - u_ex ||_{L^2}
			      loc_err[0] += std::pow(l2_error(q, u_ex, u_vec_1), 2);
			      // || u^*_k - u_ex ||_{L^2}
			      loc_err[1] += std::pow(l2_star_error(q, u_ex, u_star_1), 2);
			      // | sigma - grad(u_ex) |_{H^1}
			      loc_err[2] += std::pow(semih1_error(q, du_x_ex, du_y_ex, du_z_ex, 
			      																		u_vec_1), 2);
			      // | sigma^* - grad(u_ex) |_{H^1}
			      loc_err[3] += std::pow(semih1_star_error(q, du_x_ex,du_y_ex,du_z_ex, 
			      																					grad_star_1), 2);
			    }

			  MPI_Reduce(&local_h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);

			  MPI_Reduce(loc_err.data(), errors.data(), 4, MPI_DOUBLE, MPI_SUM,
			  						0, mpicomm);

			  for (double& d : errors)
			  	d = std::sqrt(d);

			  nnodes.push_back(tmsh.num_global_nodes());
			  nquads.push_back(num_quads);
			  h_vec.push_back(global_h);

			  MPI_Barrier(MPI_COMM_WORLD);
			  if (rank == 0)
			    {
			    	std::cout << "nnodes = " << nnodes[nnodes.size()-1] << std::endl;
			    	std::cout << "nquads = " << nquads[nquads.size()-1] << std::endl;
			    	std::cout << "h = " << h_vec[h_vec.size()-1] << std::endl;
			      for (unsigned j=0; j<4; ++j)
			        std::cout << "err" << j << " = " << errors[j] << std::endl;
			      std::cout << std::endl;
			    }

      	// ******************************************************************* 

      // Export mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }      
      sprintf(name,"p4est_refine_box_%4.4d",file_number++);
      tmsh.vtk_export (name);
      if (rank == 0) { toc ("*** Export ***"); }

    }

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }
  
  MPI_Finalize ();
  return 0;
}
