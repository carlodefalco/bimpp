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

void
lis::set_lhs_structure_aij
(int n,
 std::vector<int> &ir,
 std::vector<int> &jc)
{
  n_row = n;

  row_ptr.resize (n_row + 1);
  row_ptr.assign (n_row + 1, 0);

  jcol.resize (jc.size ());
  jcol.assign (jc.size (), 0);

  ordering_map.resize (jc.size ());
  ordering_map.assign (jc.size (), 0);
  
  for (unsigned int i = 0; i < ir.size (); ++i)
    row_ptr[ir[i] - index_base]++;

  for (int i = 0, cumsum = index_base; i < n_row; ++i)
    {
      int temp = row_ptr[i];
      row_ptr[i] = cumsum;
      cumsum += temp;
    }
  row_ptr[n_row] = ir.size () + index_base;
  for (unsigned int i = 0; i < jc.size (); ++i)
    {
      int row = ir[i];
      int dest = row_ptr[row];

      jcol[dest] = jc[i];
      ordering_map[i] = dest;
      row_ptr[row]++;
    }
  for (int i = 0, last = index_base; i <= n_row; ++i)
    {
      int temp = row_ptr[i];
      row_ptr[i] = last;
      last = temp;
    }
}

void
lis::set_lhs_structure_csr
(int n,
 std::vector<int> &ir,
 std::vector<int> &jc)
{
  n_row = n;
  row_ptr.assign (n_row + 1, 0);
  jcol.assign (jc.size (), 0);

  row_ptr = ir;
  jcol = jc;
}

int
lis::analyze_master ()
{

  map_i_s.assign (size, 0);
  map_row_s.assign (size, 0);
  map_n.assign (size, 0);
  map_nnz.assign (size, 0);

  map_n[0] = n_row / size + n_row % size;
  map_nnz[0] = row_ptr[map_n[0]] - index_base;

  for (int k = 1; k < size; ++k)
    {
      map_n[k] = n_row / size;
      map_i_s[k] = row_ptr[map_n[k] * k + n_row % size] -
        index_base;
      map_nnz[k] = row_ptr[map_n[k] * (k + 1) + n_row % size] -
        row_ptr[map_n[k] * k + n_row % size];
      map_row_s[k] = map_n[k] * k + n_row % size;

    }

  MPI_Scatter (&map_i_s[0], 1, MPI_INT,
               &i_s, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatter (&map_row_s[0], 1, MPI_INT,
               &row_s, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatter (&map_n[0], 1, MPI_INT,
               &n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Scatter (&map_nnz[0], 1, MPI_INT,
               &nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);


  MPI_Scatterv (&jcol[0], &map_nnz[0], &map_i_s[0], MPI_INT,
                MPI_IN_PLACE, 0, MPI_INT, 0, MPI_COMM_WORLD);

  for (int k = 1; k < size; ++k)
    MPI_Send (&row_ptr[map_row_s[k]], map_n[k] + 1, MPI_INT,
              k, 0, MPI_COMM_WORLD);

  return 1;
}

void
lis::set_lhs_data (std::vector<double> &xa)
{
  data = &*xa.begin ();
}


int
lis::factorize_master ()
{

  if (ordering_map.size () != 0)
    {
      double *temp = new double[ordering_map.size ()];
      for (unsigned int i = 0; i < ordering_map.size (); ++i)
        temp[i] = data[i];
      for (unsigned int i = 0; i < ordering_map.size (); ++i)
        data[ordering_map[i]] = temp[i];
      delete [] temp;
      }
  // partitioning matrix entries
  MPI_Scatterv (&data[0], &map_nnz[0], &map_i_s[0], MPI_DOUBLE,
                MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return 1;
}

int
lis::solve_master ()
{
  // Partion rhs and initial guess
  MPI_Bcast (&have_initial_guess, 1, MPI_INT, 0, MPI_COMM_WORLD);

  {
    MPI_Scatterv (&rhs[0], &map_n[0], &map_row_s[0], MPI_DOUBLE,
                  MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (have_initial_guess)
      MPI_Scatterv (&initial_guess[0], &map_n[0], &map_row_s[0],
                    MPI_DOUBLE, MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  int retval = invoke_lis_solver ();

  if (verbose)
    std::cout << std::endl
              << "Number of iterations = " << iter
              << std::endl
              << "Elapsed time = " << time << std::endl;
  
  MPI_Gatherv (MPI_IN_PLACE, 0, MPI_DOUBLE, &rhs[0],
               &map_n[0], &map_row_s[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);


  return retval;
}

void
lis::cleanup_master () { }


