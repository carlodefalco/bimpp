/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem:  
  \f[ -div (|\nabla (u)|^{p-2} \nabla (u)) = f \f]

  \f[ u = 1/q \cdot (0.5^q -
     [ (x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2)]^{q/2} \:on \:boundary \f]

  \f[ f = 3.0 \f]

  \f[ p = 3.0 \f]

  Exact Solution: 
  \f[  u = 1/q \cdot (0.5^q -
           [ (x-0.5)^2 + (y-0.5)^2 + (z-0.5)^2)]^{q/2} \f]

  Linear Solver: lis

  NonLinear Solver: adaptive_inexact_newton

  Forcing Term: forcing_type1 (0.9)
*/

#include <lis.h>
#include "mumps_class.h"
#include "lis_class.h"
#include "linear_solver.h"
#include "operators.h"
#include "mesh.h"
#include "nonlinear_solver.h"
#include "adaptive_inexact_newton_class.h"
#include "abstract_nonlinear_problem.h"
#include "abstract_forcing_term.h"
#include "forcing_class.h"
#include "bim_config.h"
#include "plaplacian_class.h"

void
run_test_problem (nonlinear_solver *solver);

int main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  linear_solver *lis_solver = new lis ();
  //linear_solver *mumps_solver = new mumps ();

  nonlinear_solver *solver =
    new adaptive_inexact_newton (lis_solver, 2);

  run_test_problem (solver);

  MPI_Finalize ();
}

void
run_test_problem (nonlinear_solver *solver)
{
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  std::vector<double> exactsolution;

  std::vector<int>    bnodes;
  std::vector<double> vnodes;

  std::vector<double> uold;

  std::vector<double> fcoeff;

  double p = 3.0;
  double q = p / (p - 1);

  abstract_nonlinear_problem *plap = new plaplacian (p);
  abstract_forcing_term *forcing = new  forcing_type1 (0.9);

  double lambda = 1.0;
  double mu = 0.5;

  if (rank == 0)
    {
      std::cout << "\n\n*****\nNon Linear Test\n*****\n";

      solver->set_output_filename ("NonLinear_Plaplacian_Test1.txt");

      std::cout << "Using solver named "
                << solver->solver_name ()
                << std::endl;

      plap->read_mesh (data_dir + "mesh_in_cube.msh");
      int nnodes = plap->msh.nnodes;

      exactsolution.resize (nnodes);
      uold.resize (nnodes);

      for (int i = 0; i < nnodes; ++i)
        {
          exactsolution[i] = (1.0 / q) * (pow (0.5, q) -
            pow ((plap->msh.p (0, i) - 0.5) * (plap->msh.p (0, i) - 0.5) +
                 (plap->msh.p (1, i) - 0.5) * (plap->msh.p (1, i) - 0.5) +
                 (plap->msh.p (2, i) - 0.5) * (plap->msh.p (2, i) - 0.5),
                   q / 2.0));

          uold[i] = exactsolution[i] * (1.0 + lambda *
                    (plap->msh.p (0, i) - mu) *
                    (plap->msh.p (1, i) - mu) *
                    (plap->msh.p (2, i) - mu));
        }

      fcoeff.resize (nnodes, 3.0);

      std::vector<int> sidelist;
      sidelist.push_back (1);
      sidelist.push_back (2);
      sidelist.push_back (3);
      sidelist.push_back (4);
      sidelist.push_back (5);
      sidelist.push_back (6);

      plap->msh.boundary_nodes (sidelist, bnodes);
      vnodes.resize (bnodes.size ());

      for (unsigned int i = 0; i < vnodes.size (); ++i)
        {
          uold[bnodes[i]] = (1.0 / q) * (pow (0.5, q) -
              pow ((plap->msh.p (0, bnodes[i]) - 0.5) *
                   (plap->msh.p (0, bnodes[i]) - 0.5) +
                   (plap->msh.p (1, bnodes[i]) - 0.5) *
                   (plap->msh.p (1, bnodes[i]) - 0.5) +
                   (plap->msh.p (2, bnodes[i]) - 0.5) *
                   (plap->msh.p (2, bnodes[i]) - 0.5),
                    q / 2.0));
          vnodes[i] = 0.0;
        }

      plap->set_exact_solution (exactsolution);
      plap->set_rhs_values (fcoeff);
      plap->set_boundary_conditions (vnodes, bnodes);
    }

  if (rank == 0)
    {
      solver->set_problem (plap);
      solver->set_forcing_term (forcing);
      solver->set_initial_guess (uold);
    }

  solver->set_max_iterations (100);
  solver->set_tolerance (1e-10);
  solver->set_min_residual (1e-10);
  solver->set_norm_type (L2);

  if (solver->linear_solver_type () == "iterative")
    {
      solver->set_max_iterations_of_linear_solver (1000);
      solver->set_iterative_method_of_linear_solver
              ("Conjugate Gradient");
      solver->set_initial_tolerance_of_linear_solver (0.5);
      solver->set_convergence_condition_of_linear_solver ("norm2_of_rhs");
    }

  bool converged = solver->solve ();

  if (rank == 0)
    {
      if (converged)
        {
          int iteration = 0;
          double residual_norm = 0.0;
          double delta_norm = 0.0;

          solver->get_result_iterations (iteration);
          solver->get_result_residual_norm (residual_norm);

          std::vector<double> delta_exact
            (uold.size ());
          for (unsigned int i = 0; i < uold.size (); ++i)
            delta_exact[i] = uold[i] - exactsolution[i];

          bim3a_norm (plap->msh, delta_exact, delta_norm, L2);
          std::cout << std::endl
                    << "Total Newton's Iterations: "
                    << iteration << std::endl
                    << "Residual Norm: " << residual_norm << std::endl
                    << "Error: " << delta_norm << std::endl;
        }
      else
        {
          std::cerr << "Not Converged!" << std::endl;
          exit (-1);
        }
    }

  solver->cleanup ();
}
