/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
     System of nonlinear equations for \f[x = (x_1, ..., x_n)^T\f] with \f[x \in \Omega \f]
  Linear Solver: lis or numps

  NonLinear Solver: backtracking_inexact_newton_example8

  Forcing Term: forcing_type3 (0.9, 2, 0.9)
*/

#include <lis.h>
#include "mumps_class.h"
#include "lis_class.h"
#include "linear_solver.h"
#include "operators.h"
#include "mesh.h"
#include "nonlinear_solver.h"
#include "backtracking_inexact_newton_example8_class.h"
#include "abstract_nonlinear_problem.h"
#include "abstract_forcing_term.h"
#include "forcing_class.h"
#include "bim_config.h"
#include "example8_class.h"
#include <sstream>
void
run_test_problem (nonlinear_solver *solver);

int main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  linear_solver *lis_solver = new lis (); // non sono riuscita a trovare il GMRES
 // linear_solver *lis_solver = new mumps ();
  
  nonlinear_solver *solver =
    new backtracking_inexact_newton_example8(lis_solver);
  
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

  std::vector<double> uold, b1,b2;

  const int n=100;
  abstract_nonlinear_problem *nonlinear_system = new example8 (n);
  
  abstract_forcing_term *forcing = new  forcing_type3 (0.9, 2, 0.9);
  
  if (rank == 0)
    {
      std::cout << "\n\n*****\nNon Linear Test\n*****\n";

      solver->set_output_filename ("Example8_Test.txt");

      std::cout << "Using solver named "
                << solver->solver_name ()
                << std::endl;


      exactsolution.assign (n , 1.0);
      uold.assign (n , 0.0);
      // Assegno i bordi del dominio
      b1.assign (n , 0.5);
      b2.assign (n , 2.0);
      b1[0] = 0.8; 
      // Assegno x0 
      for (int i=0; i <20 ; ++i)
	{
          uold[i]=0.1; 
	}
      for (int i=20; i <n ; ++i)
	{
          uold[i]=0.6;
	}
      nonlinear_system->set_exact_solution (exactsolution);
    }



  if (rank == 0)
    {
      solver->set_problem (nonlinear_system);
      solver->set_forcing_term (forcing);
      solver->set_initial_guess (uold);
      if (solver->solver_name () == "Backtracking Inexact Newton Example8")
        {
          auto tmp = static_cast<backtracking_inexact_newton_example8*> (solver);
          tmp -> set_backtracking_parameters (1e-4, 1e-4, 0, 1);
          tmp -> set_thetaPN (0.5);
          tmp -> set_thetaPG (0.8);
	  tmp -> set_backtracking_max_it(20);
	  tmp -> set_bounds (b1, b2);
      }
    }

  solver->set_max_iterations (100);
  solver->set_tolerance (1e-12);
  solver->set_min_residual (1e-12);
  solver->set_norm_type (L2);
  solver->set_max_iterations_of_linear_solver(n); 
  solver->set_initial_tolerance_of_linear_solver(.765518617913987);   
  if (solver->linear_solver_type () == "iterative")
    {
      auto tmp = static_cast<backtracking_inexact_newton_example8*> (solver);
      std::stringstream opt;
      opt << "-maxiter " << n
          << " -tol " << .765518617913987
          << " -i " << "gmres "
	  << " -restart " << n << " "
          << " -p " << "none "
          << " -conv_cond " << "norm2_r ";

      
      opt << " -initx_zeros true ";

     // std::cout << std::endl << opt.str () << std::endl;
      
      static_cast<lis*> (tmp->lin_solver)->option_string = opt.str ();
      static_cast<lis*> (tmp->lin_solver)->option_string_set = true;

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
            std::cout << std::endl
                    << "Total Newton's Iterations: "
                    << iteration << std::endl
                    << "Residual Norm: " << residual_norm << std::endl;
               //     << "Error: " << delta_norm << std::endl;
        }
      else
        {
          std::cerr << "Not Converged!" << std::endl;
          exit (-1);
        }
    }

  solver->cleanup ();
}

