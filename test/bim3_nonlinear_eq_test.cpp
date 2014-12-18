/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/
/*
  Problem:	x^2-2=0
  
  Exact Solution:	x=sqrt(2)
*/

#include <stdio.h>
#include <lis_config.h>
#include <lis.h>
#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <lis_operators.h>
#include <nonlinear_solvers.h>
#include <nonlinear_operators.h>
#include <nonlinear_examples.h>
#include <fstream>
#include <stdlib.h>
#include <cmath>

LIS_INT main(LIS_INT argc, char* argv[])
{
  int rank;
  int size;

  lis_initialize (&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
	
  std::vector<double> sol, exactsolution, uold;

  equation eq;

  inexact_newton_option option;
  option.maxIter=100;
  option.minRes=10e-10;
  option.tol=10e-10;
  option.forcing=0.5;
  option.type=Inf;

  inexact_newton_status status;
  stream_option stream = {2, "Solution_Equation.txt"};
  linear_solver_option lis_option = {"-tol 0.5", "", "-conv_cond 1 -i cg"};
  ForcingCostant forcing;

  if(rank == 0)
    {
      std::cout << "\n\n*****\nTest: Non Linear Equation\n*****\n";

      exactsolution.resize (1);
      exactsolution[0] = sqrt(2);
      uold.resize (1);
      uold[0] = 0;
      eq.import (exactsolution);
      if(stream.verbosity >= 1)
	std::cout << "\nResult of Non Linear Test \nwill be written in " << stream.filename << std::endl <<std::endl;
    }			

  status = inexact_newton<equation, ForcingCostant, LIS> (eq, uold, sol, forcing, lis_option, option, stream);
  if(rank == 0)
    {
      std::cout << "Solution: " << sol[0] << std::endl;
      std::cout << "Error: " << fabs (sol[0] - exactsolution[0]) << std::endl;
    }
  lis_finalize ();

  return 0;
}
