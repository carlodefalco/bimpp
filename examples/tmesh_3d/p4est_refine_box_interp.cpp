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


/// main
int main(int argc, char ** argv)
{
  MPI_Init (&argc, &argv);
  
  int       recursive, partforcoarsen, balance;
  MPI_Comm  mpicomm = MPI_COMM_WORLD;  							
  int       rank, size;
  tmesh_3d  tmsh;

  using q1_vec    = q1_vec<distributed_vector>;
  using idx_t     = tmesh_3d::idx_t; 

  unsigned nref_1 = 2;                // number of initial uniform refinements
  unsigned nref_2 = 4;                // number of iterative refinements

  // Mesh parameters
  std::vector<idx_t> nnodes (nref_2, 0);         // number of nodes
  std::vector<double> h_step (nref_2, 0.);       // mesh size

  // Error at every step - ||u - u_ex||_L^2(q).
  std::vector<double> error_intp (nref_2,0.);

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);
    
  MPI_Barrier (MPI_COMM_WORLD); { if (rank == 0) tic (); }
  
  // Create mesh.
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
  
  // Define function for rectangle refinement.
  std::function<int (tmesh_3d::quadrant_iterator)> box_refinement =
    [L,h] (tmesh_3d::quadrant_iterator qi)
    { return rectangle_list_refinement(qi, L, h); };

  if (rank == 0) { toc ("*** Initialization ***"); }

  // Definition of exact solution.
  func3 u_ex = 
    [L] (double x, double y, double z) -> double
    {
      return (sin(2*M_PI*x/L) * cos(2*M_PI*y/L) * sin(2*M_PI*z/L));
    };

  recursive = 0;
  partforcoarsen = 1;

  // Initial uniform refinement.
  for (int cycle = 0; cycle < nref_1; ++cycle)
    {
      // Refine mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine(recursive, partforcoarsen);        
      if (rank == 0) 
        { 
          sprintf (step, "*** Refinement and balancing %3.3d ***", 
                    cycle);
          toc (step); 
        }

      // Export refined mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      tmsh.vtk_export ((std::string("p4est_refine_box_interp_initial_mesh_")
                          + std::to_string(cycle)).c_str());
      if (rank == 0) { toc ("*** Export initial mesh ***"); }
    }
  
  // Iterative refinement.
  for (int cycle = 0; cycle < nref_2; ++cycle)
    {
      // Definition of u.
      q1_vec u_vec (tmsh.num_owned_nodes());
      bim3a_solution_with_ghosts (tmsh, u_vec, replace_op);
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

      // Export solution.
      tmsh.octbin_export((std::string("p4est_refine_box_interp_u_") 
                          + std::to_string(cycle)).c_str(), 
                          u_vec);

      // Refinement
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }

      tmsh.set_refine_marker (uniform_refinement);  // uniform refinement
      tmsh.refine (recursive, partforcoarsen);
      
      if (rank == 0) 
        { 
          sprintf (step, "*** Refinement and balancing %3.3d ***", 
                    (cycle+nref_1));
          toc (step); 
        }

      // Export refined mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      tmsh.vtk_export ((std::string("p4est_refine_box_interp_mesh_")
                          + std::to_string(cycle)).c_str());
      if (rank == 0) { toc ("*** Export ***"); }

      // Interpolate u on the new mesh
      q1_vec u_vec_intp (tmsh.num_owned_nodes());
      bim3a_solution_with_ghosts (tmsh, u_vec_intp, replace_op);
      interpolate_vector (tmsh, u_vec, u_vec_intp);

      tmsh.octbin_export((std::string("p4est_refine_box_interp_u_intp_")
                          + std::to_string(cycle)).c_str(), 
                          u_vec_intp);

      // Compute number of nodes.
      nnodes[cycle] = tmsh.num_global_nodes();

      // Compute mesh size and errors.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      double  hx = 0, hy = 0, hz = 0, h = 0;

      double err_intp = 0.0;

      for (auto q = tmsh.begin_quadrant_sweep();
          q != tmsh.end_quadrant_sweep();
          ++q)
        {
          // mesh size
          double hx = q->p(0, 7) - q->p(0, 0);
          double hy = q->p(1, 7) - q->p(1, 0);
          double hz = q->p(2, 7) - q->p(2, 0);
          h = std::max(h, std::sqrt(hx*hx + hy*hy + hz*hz));

          // ||u - u_ex||_L^2(q)
          err_intp += std::pow(l2_error(q, u_ex, u_vec_intp), 2);
        }

      // Global mesh size.
      MPI_Reduce(&h, &h_step[cycle], 1, MPI_DOUBLE, MPI_MAX, 0, mpicomm);

      // Global error - ||u - u_ex||_L^2(q).
      MPI_Reduce (&err_intp, &error_intp[cycle], 1, MPI_DOUBLE, MPI_SUM, 0, 
                  mpicomm);
      error_intp[cycle] = std::sqrt(error_intp[cycle]);

      MPI_Barrier (MPI_COMM_WORLD); 
      if (rank == 0) 
        toc ("Compute h and error ");
    }

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) 
    { 
      // timing report
      print_timing_report ();

      // mesh size, errors and estimators
      for (unsigned step = 0; step < nnodes.size(); ++step)
        {
          std::cout << "\nStep " << step << ", #nodes: "
                    << nnodes[step] << ", h: "
                    << h_step[step] << std::endl;
          std::cout << "\n\tL2 norm (intp) = " << error_intp[step] 
                    << std::endl;
          std::cout << std::endl;
        }
    }

  MPI_Finalize ();
  return 0;
}

