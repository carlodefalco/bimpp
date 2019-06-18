#include <mumps_class.h>
#include <quad_operators_3d.h>
#include <bim_distributed_vector.h>
#include <bim_sparse_distributed.h>
#include <bim_timing.h>
#include <simple_connectivity_3d.h>

#include <cassert>
#include <limits>

constexpr bool EXPORT_ALL = true;

// uniform_refinement:
//    returns 1 ----> all quadrants are refined
static int
uniform_refinement (tmesh_3d::quadrant_iterator q)
{ return 1; }

// Number of refinement steps
constexpr unsigned unif_refine_steps  = 2;  // initial uniform refinement
constexpr unsigned adapt_refine_steps = 10;  // adaptive refinement
constexpr double tol = 1.e-3;               // tolerance for refinement

// Problem parameters
constexpr double inv_epsilon = 1e11;  // 1 / epsilon
constexpr double R = 0.3;             // radius of the internal sphere
constexpr double kS = 1.;             // diffusion coefficient in the sphere
constexpr double kG = 1.;             // diffusion coefficient outside

static inline double
rho2 (double x, double y, double z)
{
  x -= .5; y -= .5; z -= .5;
  return (x*x + y*y + z*z);
}

// Exact solution
static inline double
u_ex (double x, double y, double z)
{
  double R2 = R*R;
  double r2 = rho2 (x, y, z);
  return (r2 > R2 ?  std::sin (r2) : std::sin (R2));
}

// Exact derivatives:
//
// x direction
static inline double
du_x_ex (double x, double y, double z)
{
  double R2 = R*R;
  double r2 = rho2 (x, y, z);
  return (r2 > R2 ?  2 * (x - 0.5) * std::cos (r2) : 0.);
}
//
// y direction
static inline double
du_y_ex (double x, double y, double z)
{
  double R2 = R*R;
  double r2 = rho2 (x, y, z);
  return (r2 > R2 ?  2 * (y - 0.5) * std::cos (r2) : 0.);
}
//
// z direction
static inline double
du_z_ex (double x, double y, double z)
{
  double R2 = R*R;
  double r2 = rho2 (x, y, z);
  return (r2 > R2 ?  2 * (z - 0.5) * std::cos (r2) : 0.);
}

// Load term
static inline double
load (double r2)
{
  double R2 = R*R;
  return (r2 > R2 ?
          (-6. * std::cos (r2) + 4 * r2 * std::sin (r2)) :
          inv_epsilon * std::sin (R2));
}

// Diffusion term
static inline double
diffusion (double r2)
{ return (r2 > R*R ? kG : kS); }

// Reaction term
static inline double
reaction (double r2)
{ return (r2 > R*R ? 0.0 : inv_epsilon); }

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);

  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;
  int                   rank, size;
  tmesh_3d              tmsh;

  using q1_vec          = q1_vec<distributed_vector>;
  using gradient3       = gradient3<q1_vec>;
  using idx_t           = tmesh_3d::idx_t;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  // Mesh parameters
  std::vector<idx_t> nnodes (adapt_refine_steps, 0);    // number of nodes
  std::vector<double> h_step (adapt_refine_steps, 0.);  // mesh size

  // Errors at every step
  std::vector<double> error (adapt_refine_steps,0.);  // ||u - u_ex||_L^2(q)
  std::vector<double> errorH1 (adapt_refine_steps,0.);  // |u - u_ex|_H^1(q)
  std::vector<double> errorStar (adapt_refine_steps,0.);  // ||u_star - u_ex||_L^2(q)
  std::vector<double> errorH1Star (adapt_refine_steps,0.); // ||du_star - grad(u_ex)||_L^2(q)

  // Estimators at every step
  std::vector<double> estSol (adapt_refine_steps,0.);  // ||u^* - u||_L^2(q)
  std::vector<double> estGrad (adapt_refine_steps,0.); // ||grad^* u - grad u||_L^2(q)

  // Mesh generation
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  // Initial uniform refinement
  recursive = 0; partforcoarsen = 1;
  for (unsigned cycle = 0; cycle < unif_refine_steps; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }

  // Export initial mesh
  tmsh.vtk_export ("p4est_dr_test_2_metrics_initial_mesh");

  // Adaptive refinement loop
  for (unsigned adapt = 0; adapt < adapt_refine_steps; ++adapt)
    {
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) {
        std::cout << "*** Step "  << adapt << " ***" << std::endl; }

      // Compute coefficients

      // diffusion
      std::vector<double> alpha (tmsh.num_local_quadrants (), 1.);
      q1_vec psi (tmsh.num_owned_nodes ());

      // reaction
      std::vector<double> delta (tmsh.num_local_quadrants (), 1.);
      q1_vec zeta (tmsh.num_owned_nodes ());

      // rhs
      std::vector<double> f (tmsh.num_local_quadrants (), 1.);
      q1_vec g (tmsh.num_owned_nodes ());

      double x = .0, y = .0, z = .0, r2 = .0;
      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          alpha[quadrant->get_forest_quad_idx ()] =
            diffusion (rho2 (quadrant->centroid(0),
                             quadrant->centroid(1),
                             quadrant->centroid(2)));
          for (int ii = 0; ii < 8; ++ii)
            {
              x = quadrant->p (0,ii);
              y = quadrant->p (1,ii);
              z = quadrant->p (2,ii);
              r2 = rho2 (x, y, z);
              if (! quadrant->is_hanging (ii))
                {
                  psi[quadrant->gt (ii)] = 0.;
                  zeta[quadrant->gt (ii)] = reaction (r2);
                  g[quadrant->gt (ii)] = load (r2);
                }
              else
                for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                  {
                    psi[quadrant->gparent (jj, ii)] += 0.;
                    zeta[quadrant->gparent (jj, ii)] += 0.;
                    g[quadrant->gparent (jj, ii)] += 0.;
                  }
            }
        }
      psi.assemble (replace_op);
      zeta.assemble (replace_op);
      g.assemble (replace_op);

      // Assemble system matrix and right-hand side.
      distributed_sparse_matrix A;
      A.set_ranges(tmsh.num_owned_nodes());

      // advection_diffusion
      bim3a_advection_diffusion (tmsh, alpha, psi, A);

      // reaction
      bim3a_reaction (tmsh, delta, zeta, A);

      // rhs
      q1_vec rhs (tmsh.num_owned_nodes (), mpicomm);
      bim3a_rhs (tmsh, f, g, rhs);


      // Set boundary conditions.
      dirichlet_bcs3 bcs;
      for (int i = 0; i < 6; ++i)
        bcs.push_back (std::make_tuple (0, i, u_ex));
      //
      bim3a_dirichlet_bc (tmsh, bcs, A, rhs);

      A.assemble ();
      rhs.assemble ();


      // Solve problem.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      mumps mumps_solver;
      mumps_solver.set_rhs_distributed (rhs);

      std::vector<double> vals, vals_tmp;
      std::vector<int> irow, irow_tmp, jcol, jcol_tmp;

      A.aij (vals_tmp, irow_tmp, jcol_tmp, mumps_solver.get_index_base ());
      int countzeros = std::count_if (vals_tmp.begin (), vals_tmp.end (), [] (double x) {return std::abs (x) == .0; });
      std::cout << "countzeros = " << countzeros << std::endl;
      vals.reserve (vals_tmp.size ());
      irow.reserve (irow_tmp.size ());
      jcol.reserve (jcol_tmp.size ());
      for (int ii = 0; ii < vals_tmp.size (); ++ii)
        if (std::abs (vals_tmp[ii]) > 0)
          {
            vals.push_back (vals_tmp[ii]);
            irow.push_back (irow_tmp[ii]);
            jcol.push_back (jcol_tmp[ii]);
          }
      mumps_solver.set_lhs_distributed ();
      mumps_solver.set_distributed_lhs_structure (tmsh.num_global_nodes (), irow, jcol);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("prepare solver "); }


      // Solve.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      mumps_solver.analyze ();
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("analyze "); }

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      mumps_solver.set_distributed_lhs_data (vals);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("set lhs data "); }

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      auto val = mumps_solver.factorize ();
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("factorize "); }

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      mumps_solver.solve ();
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("solve "); }


      // Get distributed solution
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      q1_vec mumps_result = mumps_solver.get_distributed_solution ();

      q1_vec result (tmsh.num_owned_nodes (), mpicomm);

      for (auto ii = result.get_range_start (); ii != result.get_range_end (); ++ii)
        result[ii] = mumps_result[ii];
      bim3a_solution_with_ghosts (tmsh, result, replace_op);
      result.assemble (replace_op);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("scatter solution "); }


      // Export solution.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_u_")
                           + std::to_string (adapt)).c_str (), result);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("export solution "); }


      // Compute exact solution on the mesh and difference between exact
      // and numerical solution
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      //
      // exact solution
      q1_vec uex (tmsh.num_owned_nodes ());
      bim3a_solution_with_ghosts (tmsh, uex);
      //
      // difference
      q1_vec diff (tmsh.num_owned_nodes ());
      bim3a_solution_with_ghosts (tmsh, diff);

      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        for (int i = 0; i < 8; ++i)
          if (! quadrant->is_hanging (i))
            {
              uex[quadrant->gt(i)] = u_ex(quadrant->p(0, i),
                                          quadrant->p(1, i),
                                          quadrant->p(2, i));
              diff[quadrant->gt(i)] = std::abs(uex[quadrant->gt(i)] - 
                                              result[quadrant->gt(i)]);
            }
      uex.assemble (replace_op);
      diff.assemble (replace_op);

      // Export exact solution and difference
      tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_uex_")
                           + std::to_string (adapt)).c_str (), uex);
      tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_diff_")
                           + std::to_string (adapt)).c_str (), diff);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) 
        { toc ("export exact solution and difference "); }


      // Cleanup
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      mumps_solver.cleanup ();
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("cleanup "); }


      // Compute reconstructed gradient.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }

      // Activation function outside the sphere
      active_fun3 regionG = [] (tmesh_3d::quadrant_iterator q)
        {
          return (rho2 (q->centroid(0),
                        q->centroid(1),
                        q->centroid(2)) > (R*R));
        };

      // Activation function inside the sphere
      active_fun3 regionS = [] (tmesh_3d::quadrant_iterator q)
        {
          return (rho2 (q->centroid(0),
                        q->centroid(1),
                        q->centroid(2)) <= (R*R));
        };

      // Gradient outside the sphere
      gradient3 du0 = bim3c_quadtree_pde_recovered_gradient
        (tmsh, result, regionG);
      // Gradient inside the sphere
      gradient3 du1 = bim3c_quadtree_pde_recovered_gradient
        (tmsh, result,  regionS);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("gradient recovery "); }

      // Reconstructed solution outside the sphere
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      q2_vec3 u_star0 = bim3c_quadtree_pde_recovered_solution
        (tmsh, result, du0);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("solution recovery outside"); }
      // Reconstructed solution inside the sphere
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      q2_vec3 u_star1 = bim3c_quadtree_pde_recovered_solution
        (tmsh, result, du1);
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("solution recovery inside"); }

      if (EXPORT_ALL)
        {
          // Export reconstructed gradients:

          MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
          // derivative along x outside the sphere
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du0_x_")
                               + std::to_string(adapt)).c_str(),
                              std::get<0>(du0));
          //
          // derivative along y outside the sphere
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du0_y_")
                               + std::to_string(adapt)).c_str(),
                              std::get<1>(du0));
          //
          // derivative along z outside the sphere
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du0_z_")
                               + std::to_string(adapt)).c_str(),
                              std::get<2>(du0));
          //
          // derivative along x inside the sphere
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du1_x_")
                               + std::to_string(adapt)).c_str(),
                              std::get<0>(du1));
          //
          // derivative along y inside the sphere
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du1_y_")
                               + std::to_string(adapt)).c_str(),
                              std::get<1>(du1));
          //
          // derivative along z inside the sphere
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du1_z_")
                               + std::to_string(adapt)).c_str(),
                              std::get<2>(du1));
          MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("export gradients "); }

          // Compute exact derivatives on the mesh and differences between 
          // exact and numerical derivatives
          MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
          //
          // derivative along x
          q1_vec duxex (tmsh.num_owned_nodes ());
          q1_vec diff_x_0 (tmsh.num_owned_nodes ());
          q1_vec diff_x_1 (tmsh.num_owned_nodes ());
          //
          // derivative along y
          q1_vec duyex (tmsh.num_owned_nodes ());
          q1_vec diff_y_0 (tmsh.num_owned_nodes ());
          q1_vec diff_y_1 (tmsh.num_owned_nodes ());
          //
          // derivative along z
          q1_vec duzex (tmsh.num_owned_nodes ());
          q1_vec diff_z_0 (tmsh.num_owned_nodes ());
          q1_vec diff_z_1 (tmsh.num_owned_nodes ());

          for (auto quadrant = tmsh.begin_quadrant_sweep ();
               quadrant != tmsh.end_quadrant_sweep ();
               ++quadrant)
            for (int i = 0; i < 8; ++i)
              if (! quadrant->is_hanging (i))
                {
                  duxex[quadrant->gt(i)] = du_x_ex(quadrant->p(0, i),
                                                    quadrant->p(1, i),
                                                    quadrant->p(2, i));
                  duyex[quadrant->gt(i)] = du_y_ex(quadrant->p(0, i),
                                                    quadrant->p(1, i),
                                                    quadrant->p(2, i));
                  duzex[quadrant->gt(i)] = du_z_ex(quadrant->p(0, i),
                                                    quadrant->p(1, i),
                                                    quadrant->p(2, i));
                  diff_x_0[quadrant->gt(i)] = std::abs(duxex[quadrant->gt(i)] - 
                                                std::get<0>(du0)[quadrant->gt(i)]);
                  diff_x_1[quadrant->gt(i)] = std::abs(duxex[quadrant->gt(i)] - 
                                                std::get<0>(du1)[quadrant->gt(i)]);
                  diff_y_0[quadrant->gt(i)] = std::abs(duyex[quadrant->gt(i)] - 
                                                std::get<1>(du0)[quadrant->gt(i)]);
                  diff_y_1[quadrant->gt(i)] = std::abs(duyex[quadrant->gt(i)] - 
                                                std::get<1>(du1)[quadrant->gt(i)]);
                  diff_z_0[quadrant->gt(i)] = std::abs(duzex[quadrant->gt(i)] - 
                                                std::get<2>(du0)[quadrant->gt(i)]);
                  diff_z_1[quadrant->gt(i)] = std::abs(duzex[quadrant->gt(i)] - 
                                                std::get<2>(du1)[quadrant->gt(i)]);
                }
              else
                {
                  duxex[quadrant->gt(i)] += 0.0;
                  duyex[quadrant->gt(i)] += 0.0;
                  duzex[quadrant->gt(i)] += 0.0;
                  diff_x_0[quadrant->gt(i)] += 0.0;
                  diff_x_1[quadrant->gt(i)] += 0.0;
                  diff_y_0[quadrant->gt(i)] += 0.0;
                  diff_y_1[quadrant->gt(i)] += 0.0;
                  diff_z_0[quadrant->gt(i)] += 0.0;
                  diff_z_1[quadrant->gt(i)] += 0.0;
                }
          duxex.assemble (replace_op);
          duyex.assemble (replace_op);
          duzex.assemble (replace_op);
          diff_x_0.assemble (replace_op);
          diff_x_1.assemble (replace_op);
          diff_y_0.assemble (replace_op);
          diff_y_1.assemble (replace_op);
          diff_z_0.assemble (replace_op);
          diff_z_1.assemble (replace_op);

          // Export exact derivatives
          //
          // derivative along x
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du_x_ex_")
                               + std::to_string(adapt)).c_str(), duxex);
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_diff0_x_")
                               + std::to_string(adapt)).c_str(), diff_x_0);
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_diff1_x_")
                               + std::to_string(adapt)).c_str(), diff_x_1);
          //
          // derivative along y
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du_y_ex_")
                               + std::to_string(adapt)).c_str(), duyex);
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_diff0_y_")
                               + std::to_string(adapt)).c_str(), diff_y_0);
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_diff1_y_")
                               + std::to_string(adapt)).c_str(), diff_y_1);
          //
          // derivative along z
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_du_z_ex_")
                               + std::to_string(adapt)).c_str(), duzex);
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_diff0_z_")
                               + std::to_string(adapt)).c_str(), diff_z_0);
          tmsh.octbin_export ((std::string ("p4est_dr_test_2_metrics_diff1_z_")
                               + std::to_string(adapt)).c_str(), diff_z_1);
          MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) 
            { toc ("export exact derivatives and differences "); }
        }

      // Solution estimator
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }

      auto estimator = [&u_star0, &u_star1, &result]
        (tmesh_3d::quadrant_iterator q)
        {
          if (rho2(q->centroid(0), q->centroid(1), q->centroid(2)) > (R*R))
            return estimator_sol (q, u_star0, result);
          else
            return estimator_sol (q, u_star1, result);
        };

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) 
        { toc ("define solution estimator "); }


      // Gradient estimator
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }

      auto grad_estimator = [&du0, &du1, &result]
        (tmesh_3d::quadrant_iterator q)
        {
          if (rho2(q->centroid(0), q->centroid(1), q->centroid(2)) > (R*R))
            return estimator_grad (q, du0, result);
          else
            return estimator_grad (q, du1, result);
        };

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) 
        { toc ("define gradient estimator "); }


      // Compute h, errors and estimators
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      double  hx = 0, hy = 0, hz = 0,
        h = std::numeric_limits<double>::max ();

      double err = 0.0;
      double errH1 = 0.0;
      double errstar = 0.0;
      double errH1star = 0.0;
      double estsol = 0.0;
      double estgrad = 0.0;

      std::vector<double> sol_est (tmsh.num_local_quadrants());
      std::vector<double> grad_est (tmsh.num_local_quadrants());

      for (auto quadrant = tmsh.begin_quadrant_sweep ();
           quadrant != tmsh.end_quadrant_sweep ();
           ++quadrant)
        {
          // ||u^* - u||_L^2(q)
          sol_est[quadrant->get_forest_quad_idx ()] = estimator(quadrant);
          estsol += std::pow(sol_est[quadrant->get_forest_quad_idx ()], 2);

          // ||grad^* u - grad u||_L^2(q)
          grad_est[quadrant->get_forest_quad_idx ()] = grad_estimator(quadrant);
          estgrad += std::pow(grad_est[quadrant->get_forest_quad_idx ()], 2);

          hx = quadrant->p(0, 1) - quadrant->p(0, 0);
          hy = quadrant->p(1, 7) - quadrant->p(1, 0);
          hy = quadrant->p(2, 7) - quadrant->p(2, 0);

          h = std::min(h, std::sqrt(hx*hx + hy*hy + hz*hz));

          // ||u - u_ex||_L^2(q)
          err += std::pow(l2_error(quadrant, u_ex, result), 2);

          // |u - u_ex|_H^1(q)
          errH1 += std::pow(semih1_error (quadrant, du_x_ex, du_y_ex, 
                                          du_z_ex, result), 2);

          if (regionG(quadrant))
            {
              // ||u_star - u_ex||_L^2(q)
              errstar += std::pow(l2_star_error(quadrant, u_ex, u_star0), 2);
              // ||du_star - grad(u_ex)||_L^2(q)
              errH1star += std::pow(semih1_star_error (quadrant, du_x_ex, 
                                                      du_y_ex, du_z_ex, 
                                                      du0), 2);
            }
          else
            {
              // ||u_star - u_ex||_L^2(q)
              errstar += std::pow(l2_star_error(quadrant, u_ex, u_star1), 2);
              // ||du_star - grad(u_ex)||_L^2(q)
              errH1star += std::pow(semih1_star_error (quadrant, du_x_ex,
                                                      du_y_ex, du_z_ex,
                                                      du1), 2);
            }
        }

      // Export estimators
      tmsh.octbin_export_quadrant ((std::string ("p4est_dr_test_2_metrics_sol_est_")
                           + std::to_string(adapt)).c_str(), sol_est);
      tmsh.octbin_export_quadrant ((std::string ("p4est_dr_test_2_metrics_grad_est_")
                           + std::to_string(adapt)).c_str(), grad_est);

      // Global mesh size
      MPI_Reduce (&h, &h_step[adapt], 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm);

      // Global errors
      //
      // ||u - u_ex||_L^2(q)
      MPI_Reduce (&err, &error[adapt], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      error[adapt] = std::sqrt(error[adapt]);
      //
      // |u - u_ex|_H^1(q)
      MPI_Reduce (&errH1, &errorH1[adapt], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      errorH1[adapt] = std::sqrt(errorH1[adapt]);
      //
      // ||u_star - u_ex||_L^2(q)
      MPI_Reduce (&errstar, &errorStar[adapt], 1, MPI_DOUBLE, MPI_SUM, 
                  0, mpicomm);
      errorStar[adapt] = std::sqrt(errorStar[adapt]);
      //
      // ||du_star - grad(u_ex)||_L^2(q)
      MPI_Reduce (&errH1star, &errorH1Star[adapt], 1, MPI_DOUBLE, MPI_SUM, 
                  0, mpicomm);
      errorH1Star[adapt] = std::sqrt(errorH1Star[adapt]);

      // Global estimators
      //
      // ||u^* - u||_L^2(q)
      MPI_Reduce (&estsol, &estSol[adapt], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      estSol[adapt] = std::sqrt(estSol[adapt]);
      //
      // ||grad^* u - grad u||_L^2(q)
      MPI_Reduce (&estgrad, &estGrad[adapt], 1, MPI_DOUBLE, MPI_SUM, 0, 
                  mpicomm);
      estGrad[adapt] = std::sqrt(estGrad[adapt]);

      // Number of mesh nodes at current step
      nnodes[adapt] = tmsh.num_global_nodes ();

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("Compute h and error "); }

      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      // Break if the number of global nodes is too large
      if (tmsh.num_global_nodes () >= 10e6)
        break;
      else if (adapt < (adapt_refine_steps - 1))
        {
          // Refine.
          tmsh.set_metrics_marker (estimator, tol*std::pow (.9, adapt), 3);
          tmsh.metrics_refine (1e3);
          std::cout << "tmsh.num_global_nodes ()= "
                    << tmsh.num_global_nodes ()
                    << std::endl;
          // Export new mesh
          tmsh.vtk_export ((std::string("p4est_dr_test_2_metrics_newmesh_")
                            + std::to_string (adapt)).c_str ());
          if (tmsh.num_global_nodes () >= 10e6)
            {
              std::cout << "too many nodes!" << std::endl;
              break;
            }
        }
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc ("refine "); }
    }

  if (rank == 0)
    for (unsigned step = 0; step < nnodes.size(); ++step)
      {
        std::cout << "Step " << step << ", #nodes: "
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

  MPI_Finalize ();
  return 0;

}
