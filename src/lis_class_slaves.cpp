/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file lis_class.cpp
  \brief interface for linear solver built with for lis library.
*/

#include "lis_class.h"
#include <stdlib.h>
#include <sstream>
#include <string>
#include <cstring>


int
lis::analyze_slaves ()
{
  MPI_Scatter (&map_i_s[0], 1, MPI_INT,
               &i_s, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatter (&map_row_s[0], 1, MPI_INT,
               &row_s, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatter (&map_n[0], 1, MPI_INT,
               &n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatter (&map_nnz[0], 1, MPI_INT,
               &nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  jcol.assign (nnz, 0);
  MPI_Scatterv (&jcol[0], &map_nnz[0], &map_i_s[0], MPI_INT,
                &jcol[0], nnz, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Status *status = NULL;
  row_ptr.assign (n + 1, 0);
  MPI_Recv (&row_ptr[0], n + 1, MPI_INT,
            0, 0, MPI_COMM_WORLD, status);

  data = new double[nnz];
  rhs = new double[n];

  return 1;
}


int
lis::factorize_slaves ()
{
  // partitioning matrix entries
  MPI_Scatterv (&data[0], &map_nnz[0], &map_i_s[0], MPI_DOUBLE,
                &data[0], nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return 1;
}

int
lis::solve_slaves ()
{
  // Partion rhs and initial guess
  MPI_Bcast (&have_initial_guess, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatterv (&rhs[0], &map_n[0], &map_row_s[0], MPI_DOUBLE,
                &rhs[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (have_initial_guess)
    {
      initial_guess = new double[n];
      
      MPI_Scatterv (&initial_guess[0], &map_n[0], &map_row_s[0],
                    MPI_DOUBLE, &initial_guess[0], n, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
    }

  int retval = invoke_lis_solver ();

  MPI_Gatherv (&rhs[0], n, MPI_DOUBLE, &rhs[0],
               &map_n[0], &map_row_s[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return retval;
}

void
lis::cleanup_slaves ()
{
  delete [] data;
  delete [] rhs;
 
  if (have_initial_guess)
    delete [] initial_guess;
}
