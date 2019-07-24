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
  using gradient3 = gradient3<distributed_vector>;
  using idx_t     = tmesh_3d::idx_t;

  unsigned nref_1 = 1;                // number of initial uniform refinements
  unsigned nref_2 = 5;                // number of iterative refinements

  // Mesh parameters
  std::vector<idx_t> nnodes (nref_2, 0);      // number of nodes
  std::vector<double> h_step (nref_2, 0.);    // mesh size

  // Errors at every step
  std::vector<double> error (nref_2,0.);      // ||u - u_ex||_L^2(q)
  std::vector<double> errorH1 (nref_2,0.);    // |u - u_ex|_H^1(q)
  std::vector<double> errorStar (nref_2,0.);  // ||u_star - u_ex||_L^2(q)
  std::vector<double> errorH1Star(nref_2,0.); // ||du_star - grad(u_ex)||_L^2(q)

  // Estimators at every step
  std::vector<double> estSol (nref_2,0.);     // ||u^* - u||_L^2(q)
  std::vector<double> estGrad (nref_2,0.);    // ||grad^* u - grad u||_L^2(q)

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

  // Definition of exact derivatives.
  func3 du_x_ex =
    [L] (double x, double y, double z) -> double
    {
      return (2*M_PI *cos(2*M_PI*x/L) *cos(2*M_PI*y/L) *sin(2*M_PI*z/L) /L);
    };

  func3 du_y_ex =
    [L] (double x, double y, double z) -> double
    {
      return (-2*M_PI *sin(2*M_PI*x/L) *sin(2*M_PI*y/L) *sin(2*M_PI*z/L) /L);
    };

  func3 du_z_ex =
    [L] (double x, double y, double z) -> double
    {
      return (2*M_PI *sin(2*M_PI*x/L) *cos(2*M_PI*y/L) *cos(2*M_PI*z/L) /L);
    };

  recursive = 0;
  partforcoarsen = 1;

  // Initial uniform refinement.
  for (int cycle = 0; cycle < nref_1; ++cycle)
    {
      // Refine mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
      if (rank == 0)
        {
          sprintf (step, "*** Refinement and balancing %3.3d ***",
                   cycle);
          toc (step);
        }

      // Export refined mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      tmsh.vtk_export ((std::string("p4est_refine_box_inital_mesh_")
                        + std::to_string(cycle)).c_str());
      if (rank == 0) { toc ("*** Export initial mesh ***"); }
    }

  // Iterative refinement according to boxes.
  for (int cycle = 0; cycle < nref_2; ++cycle)
    {
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }

      tmsh.set_refine_marker (uniform_refinement);  // uniform refinement
      tmsh.refine (recursive, partforcoarsen);

      tmsh.set_refine_marker (box_refinement);      // rectangle refinement
      tmsh.refine (recursive, partforcoarsen);

      if (rank == 0)
        {
          sprintf (step, "*** Refinement and balancing %3.3d ***",
                   (cycle+nref_1));
          toc (step);
        }

      // Export refined mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      tmsh.vtk_export ((std::string("p4est_refine_box_mesh_")
                        + std::to_string(cycle)).c_str());
      if (rank == 0) { toc ("*** Export ***"); }

      // Compute number of nodes.
      nnodes[cycle] = tmsh.num_global_nodes ();

      // Definition of u.
      q1_vec u_vec (tmsh.num_owned_nodes());
      bim3a_solution_with_ghosts (tmsh,u_vec);
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

      // Compute recovered gradient.
      gradient3 grad_star = bim3c_quadtree_pde_recovered_gradient (tmsh,
                                                                   u_vec);

      // Compute recovered solution.
      q2_vec3 u_star = bim3c_quadtree_pde_recovered_solution(tmsh,u_vec,
                                                             grad_star);

      // Export solution.
      tmsh.octbin_export((std::string("p4est_refine_box_u_")
                          + std::to_string(cycle)).c_str(),
                         u_vec);

      // Compute mesh size, errors and estimators.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }

      double  hx = 0, hy = 0, hz = 0, h = 0;
      double err = 0.0;
      double errH1 = 0.0;
      double errstar = 0.0;
      double errH1star = 0.0;
      double estsol = 0.0;
      double estgrad = 0.0;

      std::vector<double> sol_est (tmsh.num_local_quadrants());
      std::vector<double> grad_est (tmsh.num_local_quadrants());

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
          err += std::pow(l2_error(q, u_ex, u_vec), 2);

          // ||u_star - u_ex||_L^2(q)
          errstar += std::pow(l2_star_error(q, u_ex, u_star), 2);

          // |u - u_ex|_H^1(q)
          errH1 += std::pow(semih1_error(q, du_x_ex, du_y_ex, du_z_ex,
                                         u_vec), 2);

          // ||du_star - grad(u_ex)||_L^2(q)
          errH1star = std::pow(semih1_star_error(q, du_x_ex,du_y_ex,du_z_ex,
                                                 grad_star), 2);

          double temp_est = 0.;

          // ||u^* - u||_L^2(q)
          temp_est = std::pow(estimator_sol(q, u_star, u_vec), 2);
          estsol += temp_est;
          sol_est[q->get_forest_quad_idx ()] = temp_est;

          // ||grad^* u - grad u||_L^2(q)
          temp_est = std::pow(estimator_grad(q, grad_star, u_vec), 2);
          estgrad += temp_est;
          grad_est[q->get_forest_quad_idx ()] = temp_est;

        }
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("Compute h and error"); }

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      // Export estimators.
      tmsh.octbin_export_quadrant ((std::string("p4est_refine_box_err_")
                                    + std::to_string(cycle)).c_str(),
                                   sol_est);
      tmsh.octbin_export_quadrant ((std::string("p4est_refine_box_err_")
                                    + std::to_string(cycle)).c_str(),
                                   grad_est);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("export error plots"); }

      // Global mesh size.
      MPI_Reduce(&h, &h_step[cycle], 1, MPI_DOUBLE, MPI_MAX, 0, mpicomm);

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      // Global errors
      //
      // ||u - u_ex||_L^2(q)
      MPI_Reduce (&err, &error[cycle], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      error[cycle] = std::sqrt(error[cycle]);
      //
      // |u - u_ex|_H^1(q)
      MPI_Reduce (&errH1, &errorH1[cycle], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      errorH1[cycle] = std::sqrt(errorH1[cycle]);
      //
      // ||u_star - u_ex||_L^2(q)
      MPI_Reduce (&errstar, &errorStar[cycle], 1, MPI_DOUBLE, MPI_SUM,
                  0, mpicomm);
      errorStar[cycle] = std::sqrt(errorStar[cycle]);
      //
      // ||du_star - grad(u_ex)||_L^2(q)
      MPI_Reduce (&errH1star, &errorH1Star[cycle], 1, MPI_DOUBLE, MPI_SUM,
                  0, mpicomm);
      errorH1Star[cycle] = std::sqrt(errorH1Star[cycle]);

      // Global estimators
      //
      // ||u^* - u||_L^2(q)
      MPI_Reduce (&estsol, &estSol[cycle], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      estSol[cycle] = std::sqrt(estSol[cycle]);
      //
      // ||grad^* u - grad u||_L^2(q)
      MPI_Reduce (&estgrad, &estGrad[cycle], 1, MPI_DOUBLE, MPI_SUM, 0,
                  mpicomm);
      estGrad[cycle] = std::sqrt(estGrad[cycle]);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("reduce global values"); }

    }

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    {
      // timing report
      print_timing_report ();

      // mesh size, errors and estimators
      for (unsigned step = 0; step < nnodes.size (); ++step)
        {
          std::cout << "\nStep " << step << ", #nodes: "
                    << nnodes[step] << ", h: "
                    << h_step[step] << std::endl;
          std::cout << "\tL2 norm = " << error[step]
                    << "\n\tH1 seminorm = " << errorH1[step]
                    << "\n\tL2* norm = " << errorStar[step]
                    << "\n\tL2* norm (gradient) = " << errorH1Star[step]
                    << "\n\tSolution estimator = " << estSol[step]
                    << "\n\tGradient estimator = " << estGrad[step]
                    << std::endl;
          std::cout << std::endl;
        }
    }

  MPI_Finalize ();
  return 0;
}

