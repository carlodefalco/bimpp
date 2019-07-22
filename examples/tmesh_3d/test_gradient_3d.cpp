#include <bim_timing.h>
#include <tmesh_3d.h>

#include <simple_connectivity_3d.h>   // one tree

#include <quad_operators_3d.h>

static constexpr unsigned nref_uniform =  2;
static constexpr unsigned nref_firstquad =  3;

// uniform_refinement:
//
// returns 1 ----> all quadrants are refined
static int
uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
{ return 1; }

// firstquad_refinement:
//
// returns 1 if the global index of node 0 is equal to 0
// ----> only the first quadrant is refined
static int
firstquad_refinement (tmesh_3d::quadrant_iterator quadrant)
{
  bool is_zero = false;
  if (quadrant->gt(0) == 0) 
    is_zero = true;
  return is_zero;
}

// f_uex:
//
// implements the function u(x,y,z)
static inline double
f_uex (double x, double y, double z)
{
  return (x + y + z +
          x*y + x*z + y*z +
          x*y*z);
}

// f_dudx_ex:
//
// implements the derivative along x dudx(x,y,z)
static inline double
f_dudx_ex (double x, double y, double z)
{
  return (1 + y + z + y*z);
}

// f_dudy_ex:
//
// implements the derivative along y dudy(x,y,z)
static inline double
f_dudy_ex (double x, double y, double z)
{
  return (1 + x + z + x*z);
}

// f_dudz_ex:
//
// implements the derivative along z dudz(x,y,z)
static inline double
f_dudz_ex (double x, double y, double z)
{
  return (1 + x + y + x*y);
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

  using q1_vec = q1_vec<distributed_vector>;
  using gradient3 = gradient3<distributed_vector>;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  // Errors at every step
  //
  // ||u - f_uex||_L^2(q)
  std::vector<double> error (nref_uniform+nref_firstquad,0.);
  // |u - f_uex|_H^1(q)
  std::vector<double> errorH1 (nref_uniform+nref_firstquad,0.);
  // ||u_star - f_uex||_L^2(q)
  std::vector<double> errorStar (nref_uniform+nref_firstquad,0.);
  // ||du_star - grad(f_uex)||_L^2(q)
  std::vector<double> errorH1Star (nref_uniform+nref_firstquad,0.);

  // Estimators at every step
  //
  // ||u^* - u||_L^2(q)
  std::vector<double> estSol (nref_uniform+nref_firstquad,0.);
  // ||grad^* u - grad u||_L^2(q)
  std::vector<double> estGrad (nref_uniform+nref_firstquad,0.);

  // Maximum differences at every step
  std::vector<double> diff_u (nref_uniform+nref_firstquad,0.);
  std::vector<double> diff_dudx (nref_uniform+nref_firstquad,0.);
  std::vector<double> diff_dudy (nref_uniform+nref_firstquad,0.);
  std::vector<double> diff_dudz (nref_uniform+nref_firstquad,0.);

  /************************ initialization of the mesh ***********************/

  MPI_Barrier (mpicomm);
  if (rank == 0) { tic (); }
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  if (rank == 0) { toc ("read connectivity"); }
  MPI_Barrier (mpicomm);

  /**************************** uniform refinement ***************************/
 
  for (unsigned iter = 0; iter < nref_uniform; ++iter)
    {
      MPI_Barrier (mpicomm);
      if (rank == 0) { tic (); }

      tmsh.set_refine_marker (uniform_refinement);
      recursive = 0; partforcoarsen = 1;
      tmsh.refine (recursive, partforcoarsen);

      // Export refined mesh.
      tmsh.vtk_export ((std::string("test_gradient_3d_uniform_refinement_")
                          + std::to_string(iter)).c_str());

      if (rank == 0) { toc ("uniform refinement"); }
      MPI_Barrier (mpicomm);

      // update of u and exact derivatives (uniform refinement)
      q1_vec uex_ur(tmsh.num_owned_nodes());
      bim3a_solution_with_ghosts(tmsh,uex_ur);

      for (auto q = tmsh.begin_quadrant_sweep();
            q != tmsh.end_quadrant_sweep();
            ++q)
        {
          for (int nn = 0; nn < 8; ++nn)
            {
              // assemble non-hanging nodes
              if (! q->is_hanging(nn))
                {
                  double x = q->p(0,nn);
                  double y = q->p(1,nn);
                  double z = q->p(2,nn);
                  uex_ur[q->gt(nn)] = f_uex(x,y,z);
                }
            }
        }
      uex_ur.assemble(replace_op);

      // computation of the recovered gradient (uniform refinement)
      gradient3 grad_ur = bim3c_quadtree_pde_recovered_gradient (tmsh,uex_ur);
      q1_vec dudx_ur = std::get<0>(grad_ur);
      q1_vec dudy_ur = std::get<1>(grad_ur);
      q1_vec dudz_ur = std::get<2>(grad_ur);

      // computation of the recovered solution (uniform refinement)
      q2_vec3 ustar_ur = bim3c_quadtree_pde_recovered_solution (tmsh,uex_ur,
                                                                grad_ur);

      // computation of errors, estimators and differences 
      // (firstquad_refinement)
      double err = 0.0, errH1 = 0.0, errstar = 0.0, errH1star = 0.0;
      double estgrad = 0.0, estsol = 0.0;
      double diffu = 0.0, diffdudx = 0.0, diffdudy = 0.0, diffdudz = 0.0;
      for (auto q = tmsh.begin_quadrant_sweep();
                q != tmsh.end_quadrant_sweep();
                ++q)
        {
          // Errors
          //
          // ||u - f_uex||_L^2(q)
          err += std::pow(l2_error(q, f_uex, uex_ur), 2);
          // |u - f_uex|_H^1(q)
          errH1 += std::pow(semih1_error (q, f_dudx_ex, f_dudy_ex, 
                                          f_dudz_ex, uex_ur), 2);
          // ||u_star - f_uex||_L^2(q)
          errstar += std::pow(l2_star_error(q, f_uex, ustar_ur), 2);
          // ||du_star - grad(f_uex)||_L^2(q)
          errH1star += std::pow(semih1_star_error (q, f_dudx_ex, 
                                                   f_dudy_ex, f_dudz_ex,
                                                   grad_ur), 2);

          // Estimators
          //
          // ||grad^* u - grad u||_L^2(q)
          estgrad += std::pow(estimator_grad(q, grad_ur, uex_ur),2);
          // ||u^* - u||_L^2(q)
          estsol += std::pow(estimator_sol(q, ustar_ur, uex_ur),2);

          // Differences
          //
          for (int nn = 0; nn < 8; ++nn)
            {
              if (! q->is_hanging(nn))
                {
                  double x = q->p(0,nn);
                  double y = q->p(1,nn);
                  double z = q->p(2,nn);

                  double dudx = f_dudx_ex(x,y,z);
                  double dudy = f_dudy_ex(x,y,z);
                  double dudz = f_dudz_ex(x,y,z);

                  diffu = std::max(std::abs(ustar_ur[q->get_forest_quad_idx ()]
                                                    [q->gt(nn)] - 
                                            uex_ur[q->gt(nn)]),
                                    diffu);

                  diffdudx = std::max(std::abs(dudx - dudx_ur[q->gt(nn)]), 
                                      diffdudx);
                  diffdudy = std::max(std::abs(dudy - dudy_ur[q->gt(nn)]), 
                                      diffdudy);
                  diffdudz = std::max(std::abs(dudz - dudz_ur[q->gt(nn)]), 
                                      diffdudz);
                }
            }
        }

      // Global errors
      //
      // ||u - f_uex||_L^2(q)
      MPI_Reduce (&err, &error[iter], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      error[iter] = std::sqrt(error[iter]);
      //
      // |u - f_uex|_H^1(q)
      MPI_Reduce (&errH1, &errorH1[iter], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      errorH1[iter] = std::sqrt(errorH1[iter]);
      //
      // ||u_star - f_uex||_L^2(q)
      MPI_Reduce (&errstar, &errorStar[iter], 1, MPI_DOUBLE, MPI_SUM, 
                  0, mpicomm);
      errorStar[iter] = std::sqrt(errorStar[iter]);
      //
      // ||du_star - grad(f_uex)||_L^2(q)
      MPI_Reduce (&errH1star, &errorH1Star[iter], 1, MPI_DOUBLE, MPI_SUM, 
                  0, mpicomm);
      errorH1Star[iter] = std::sqrt(errorH1Star[iter]);

      // Global estimators
      //
      // ||u^* - u||_L^2(q)
      MPI_Reduce (&estsol, &estSol[iter], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      estSol[iter] = std::sqrt(estSol[iter]);
      //
      // ||grad^* u - grad u||_L^2(q)
      MPI_Reduce (&estgrad, &estGrad[iter], 1, MPI_DOUBLE, MPI_SUM, 0, 
                  mpicomm);
      estGrad[iter] = std::sqrt(estGrad[iter]);

      // Global differences
      MPI_Reduce (&diffdudx, &diff_dudx[iter], 1, MPI_DOUBLE, MPI_MAX, 0, 
                  mpicomm);
      //
      MPI_Reduce (&diffdudy, &diff_dudy[iter], 1, MPI_DOUBLE, MPI_MAX, 0, 
                  mpicomm);
      //
      MPI_Reduce (&diffdudz, &diff_dudz[iter], 1, MPI_DOUBLE, MPI_MAX, 0, 
                  mpicomm);
    }


  /*************************** firstquad refinement **************************/

  for (unsigned iter = 0; iter < nref_firstquad; ++iter)
    {
      MPI_Barrier (mpicomm);
      if (rank == 0) { tic (); }

      tmsh.set_refine_marker (firstquad_refinement);
      recursive = 1; partforcoarsen = 1;
      tmsh.refine (recursive, partforcoarsen);

      // Export refined mesh.
      tmsh.vtk_export ((std::string("test_gradient_3d_firstquad_refinement_")
                          + std::to_string(iter)).c_str());

      if (rank == 0) { toc ("firstquad refinement"); }
      MPI_Barrier (mpicomm);

      // update of u (firstquad_refinement)
      q1_vec uex_fqr(tmsh.num_owned_nodes());
      bim3a_solution_with_ghosts(tmsh,uex_fqr);
      for (auto q = tmsh.begin_quadrant_sweep();
                q != tmsh.end_quadrant_sweep();
                ++q)
        {
          for (int nn = 0; nn < 8; ++nn)
            {
              // assemble non-hanging nodes
              if (! q->is_hanging(nn))
                {
                  double x = q->p(0,nn);
                  double y = q->p(1,nn);
                  double z = q->p(2,nn);
                  uex_fqr[q->gt(nn)] = f_uex(x,y,z);
                }
            }  
        }
      uex_fqr.assemble(replace_op);     

      // computation of the recovered gradient (firstquad_refinement)
      gradient3 grad_fqr = bim3c_quadtree_pde_recovered_gradient(tmsh,uex_fqr);
      q1_vec dudx_fqr = std::get<0>(grad_fqr);
      q1_vec dudy_fqr = std::get<1>(grad_fqr);
      q1_vec dudz_fqr = std::get<2>(grad_fqr);

      // computation of the recovered solution (firstquad_refinement)
      q2_vec3 ustar_fqr = bim3c_quadtree_pde_recovered_solution (tmsh,uex_fqr,
                                                                 grad_fqr);

      // computation of errors, estimators and differences 
      // (firstquad_refinement)
      double err = 0.0, errH1 = 0.0, errstar = 0.0, errH1star = 0.0;
      double estgrad = 0.0, estsol = 0.0;
      double diffu = 0.0, diffdudx = 0.0, diffdudy = 0.0, diffdudz = 0.0;
      for (auto q = tmsh.begin_quadrant_sweep();
                q != tmsh.end_quadrant_sweep();
                ++q)
        {
          // Errors
          //
          // ||u - f_uex||_L^2(q)
          err += std::pow(l2_error(q, f_uex, uex_fqr), 2);
          // |u - f_uex|_H^1(q)
          errH1 += std::pow(semih1_error (q, f_dudx_ex, f_dudy_ex, 
                                          f_dudz_ex, uex_fqr), 2);
          // ||u_star - f_uex||_L^2(q)
          errstar += std::pow(l2_star_error(q, f_uex, ustar_fqr), 2);
          // ||du_star - grad(f_uex)||_L^2(q)
          errH1star += std::pow(semih1_star_error (q, f_dudx_ex, 
                                                   f_dudy_ex, f_dudz_ex,
                                                   grad_fqr), 2);

          // Estimators
          //
          // ||grad^* u - grad u||_L^2(q)
          estgrad += std::pow(estimator_grad(q, grad_fqr, uex_fqr),2);
          // ||u^* - u||_L^2(q)
          estsol += std::pow(estimator_sol(q, ustar_fqr, uex_fqr),2);

          // Differences
          //
          for (int nn = 0; nn < 8; ++nn)
            {
              if (! q->is_hanging(nn))
                {
                  double x = q->p(0,nn);
                  double y = q->p(1,nn);
                  double z = q->p(2,nn);

                  double dudx = f_dudx_ex(x,y,z);
                  double dudy = f_dudy_ex(x,y,z);
                  double dudz = f_dudz_ex(x,y,z);

                  diffu = std::max(std::abs(ustar_fqr[q->get_forest_quad_idx()]
                                                     [q->gt(nn)] - 
                                            uex_fqr[q->gt(nn)]),
                                    diffu);

                  diffdudx = std::max(std::abs(dudx - dudx_fqr[q->gt(nn)]), 
                                      diffdudx);
                  diffdudy = std::max(std::abs(dudy - dudy_fqr[q->gt(nn)]), 
                                      diffdudy);
                  diffdudz = std::max(std::abs(dudz - dudz_fqr[q->gt(nn)]), 
                                      diffdudz);
                }
            }
        }

      unsigned iter2 = iter + nref_uniform;

      // Global errors
      //
      // ||u - f_uex||_L^2(q)
      MPI_Reduce (&err, &error[iter2], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      error[iter2] = std::sqrt(error[iter2]);
      //
      // |u - f_uex|_H^1(q)
      MPI_Reduce (&errH1, &errorH1[iter2], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      errorH1[iter2] = std::sqrt(errorH1[iter2]);
      //
      // ||u_star - f_uex||_L^2(q)
      MPI_Reduce (&errstar, &errorStar[iter2], 1, MPI_DOUBLE, MPI_SUM, 0, 
                  mpicomm);
      errorStar[iter2] = std::sqrt(errorStar[iter2]);
      //
      // ||du_star - grad(f_uex)||_L^2(q)
      MPI_Reduce (&errH1star, &errorH1Star[iter2], 1, MPI_DOUBLE, MPI_SUM, 0, 
                  mpicomm);
      errorH1Star[iter2] = std::sqrt(errorH1Star[iter2]);

      // Global estimators
      //
      // ||u^* - u||_L^2(q)
      MPI_Reduce (&estsol, &estSol[iter2], 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm);
      estSol[iter2] = std::sqrt(estSol[iter2]);
      //
      // ||grad^* u - grad u||_L^2(q)
      MPI_Reduce (&estgrad, &estGrad[iter2], 1, MPI_DOUBLE, MPI_SUM, 0, 
                  mpicomm);
      estGrad[iter2] = std::sqrt(estGrad[iter2]);

      // Global differences
      MPI_Reduce (&diffdudx, &diff_dudx[iter2], 1, MPI_DOUBLE, MPI_MAX, 0, 
                  mpicomm);
      //
      MPI_Reduce (&diffdudy, &diff_dudy[iter2], 1, MPI_DOUBLE, MPI_MAX, 0, 
                  mpicomm);
      //
      MPI_Reduce (&diffdudz, &diff_dudz[iter2], 1, MPI_DOUBLE, MPI_MAX, 0, 
                  mpicomm);
    }

  if (rank == 0) {print_timing_report();}
  MPI_Barrier (mpicomm);

  /************************* print errors & estimators ***********************/

  if (rank == 0)
    {
      std::cout << "### Uniform Refinement" << std::endl;
      for (unsigned step = 0; step < error.size(); ++step)
        {
          if (step == nref_uniform)
            std::cout << "### Firstquad Refinement" << std::endl;
          std::cout << "Step " << step << std::endl;
          std::cout << "\tL2 norm = " << error[step] 
                    << "\n\tH1 seminorm = " << errorH1[step] 
                    << "\n\tL2* norm = " << errorStar[step] 
                    << "\n\tL2* norm (gradient) = " << errorH1Star[step]
                    << "\n\tSolution estimator = " << estSol[step]
                    << "\n\tGradient estimator = " << estGrad[step]
                    << "\n\tMaximum differences: "
                    << "\n\t\tu: " << diff_u[step]
                    << "\n\t\tdudx: " << diff_dudx[step]
                    << "\n\t\tdudy: " << diff_dudy[step]
                    << "\n\t\tdudz: " << diff_dudz[step]
                    << std::endl;
          std::cout << std::endl;
        }
    }
  
  MPI_Finalize ();
  return 0;

}