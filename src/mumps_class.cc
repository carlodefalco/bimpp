/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms 
  the terms of the GNU/GPL licence v3
*/


#include <mumps_class.h>

void mumps::init ()
{
  id.job =  JOB_INIT;
  id.par =   1;                      // host working
  id.sym =   0;                      // non symmetric
  id.comm_fortran = F77_COMM_WORLD;  // MPI_COMM_WORLD
  
  dmumps_c (&id);
  
  // streams    
  id.icntl[0] =   -1; // Output stream for error messages
  id.icntl[1] =   -1; // Output stream for diagnostic messages
  id.icntl[2] =   -1; // Output stream for global information
  id.icntl[3] =    0; // Level of printing 

  // Matrix input format
  id.icntl[4]  =   0; 
  id.icntl[17] =   0; 

  // ordering
  id.icntl[6] =  0; // metis (5), or pord (4), or AMD (0), AMF (2), QAMD (6)	  
}


void mumps::set_lhs_structure (int n, std::vector<int> &ir, std::vector<int> &jc)
{
  id.n  = n;
  id.nz = ir.size ();
  
  id.irn = &*ir.begin ();
  id.jcn = &*jc.begin ();
}

void mumps::analyze ()
{
  id.job = JOB_ANALYZE;
  dmumps_c (&id);
}

void mumps::set_lhs_data (std::vector<double> &xa)
{
  // Define LHS entries
  id.a   = &*xa.begin ();
}

void mumps::set_rhs (std::vector<double> &rhs)
{
  // Define RHS 
  id.rhs  =  &*rhs.begin ();
  id.nrhs =  1;
  id.lrhs =  rhs.size ();
}

void mumps::factorize ()
{
  id.job = JOB_FACTORIZE;
  dmumps_c (&id);
}

void mumps::solve ()
{
  id.job = JOB_SOLVE;
  dmumps_c (&id);
}

void mumps::cleanup ()
{
  // clean up
  id.job =  JOB_END; 
  dmumps_c (&id);
}



