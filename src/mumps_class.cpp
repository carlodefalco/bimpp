/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file mumps_class.cpp
  \brief wrapper for mumps data.
*/

#include "mumps_class.h"

void
mumps::init ()
{
  id.job =  JOB_INIT;
  id.par =   working_host;           // host working
  id.sym =   0;                      // non symmetric
  id.comm_fortran = F77_COMM_WORLD;  // MPI_COMM_WORLD

  dmumps_c (&id);

  // streams
  int mult = verbose ? 1 : -1;
  id.icntl[0] = mult * 6; // Output stream for error messages
  id.icntl[1] = mult * 6; // Output stream for diagnostic messages
  id.icntl[2] = mult * 6; // Output stream for global information
  id.icntl[3] = mult * 3; // Level of printing

  // Matrix input format
  id.icntl[4]  =   0;
  id.icntl[17] =   0;

  // ordering
  id.icntl[6] =  7; // scotch (3), or pord (4), or metis (5), or AMD (0), AMF (2), QAMD (6), AUTO (7)

  // space for fill-in ---
  id.icntl[13] = 300;
  id.icntl[22] = (icntl23 > 0) ? icntl23 : 0;

  // iterative refinement
  id.icntl[9] = 25;
}


void
mumps::set_lhs_structure
(int n,
 std::vector<int> &ir,
 std::vector<int> &jc,
 matrix_format_t f)
{
  if (f == csr)
    {
      id.n  = n;
      id.nz = jc.size ();
      id.irn = new int[jc.size ()];
      id.jcn = &*jc.begin ();
      for (int i = 0; i < n; ++i)
        {
          for (int j = ir[i]; j < ir[i+1]; ++j)
            id.irn[j - index_base] = i + index_base;
        }
    }
  else
    {
      id.n  = n;
      id.nz = ir.size ();
      id.irn = &*ir.begin ();
      id.jcn = &*jc.begin ();
    }
}

int
mumps::analyze ()
{
  id.job = JOB_ANALYZE;
  dmumps_c (&id);
  return id.info[0];
}

void
mumps::set_lhs_data (std::vector<double> &xa)
{
  // Define LHS entries
  id.a   = &*xa.begin ();
}

void
mumps::set_rhs (std::vector<double> &rhs)
{
  // Define RHS
  id.rhs  =  &*rhs.begin ();
  id.nrhs =  1;
  id.lrhs =  rhs.size ();
}

int
mumps::factorize ()
{
  id.job = JOB_FACTORIZE;
  dmumps_c (&id);
  return id.info[0];
}

int
mumps::solve ()
{
  id.job = JOB_SOLVE;
  dmumps_c (&id);
  return id.info[0];
}

void
mumps::cleanup ()
{
  // clean up
  id.job =  JOB_END;
  dmumps_c (&id);
}
