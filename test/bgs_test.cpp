/*
  Copyright (C) 2015 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <bim_sparse.h>
#include <bim_timing.h>
#include <bgs_class.h>
#include <mumps_class.h>
#include <lis_class.h>
#include <mpi.h>
#include <fstream>

int main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  std::vector<int>   i,j;
  std::vector<double>  a;
  std::vector<double>  rhs;
  sparse_matrix sp;
  
  if (rank == 0)
    {
      tic ();
      sp.resize (21000);
      
      for (unsigned int ii = 0; ii < sp.rows (); ++ii)
        {
          unsigned int jstart, jend;
          unsigned int bw = 3;
          if (ii - bw >= 0) 
            jstart = ii - bw; 
          else 
            jstart = 0;
          if ((ii + bw + 1) <= sp.rows ()) 
            jend = ii + bw; 
          else 
            jend = sp.rows () - 1;
          for (unsigned int jj = jstart; jj <= jend; ++jj)
            if (ii == jj)
              sp[ii][jj] = 7.0;
            else
              sp[ii][jj] = -1.0;
        }
      
      sp.aij (a, i, j);
      std::cout << "rank = " << rank << std::endl;
      toc ("build matrix");
      rhs = std::vector<double> (sp.rows (), 1.0);
    }

  MPI_Barrier (MPI_COMM_WORLD);

  std::vector<linear_solver*> block_solvers;
  block_solvers.push_back (new mumps);
  block_solvers.push_back (new lis);
  ((lis *)block_solvers[1])->verbose = false;
  block_solvers.push_back (new lis);
  ((lis *)block_solvers[2])->verbose = false;
  
  bgs bgs_solver (block_solvers);

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    {
      tic ();
      bgs_solver.set_lhs_structure (sp.rows (), i, j, linear_solver::aij);
      std::cout << "rank = " << rank << std::endl;
      toc ("set lhs struct");
    }
  MPI_Barrier (MPI_COMM_WORLD);
  
  if (rank == 0)
    tic ();

  bgs_solver.analyze ();

  if (rank == 0)
    {
      std::cout << "rank = " << rank << std::endl;
      toc ("analyze");
    }
  
  // first solve
  if (rank == 0)
    {
      tic ();
      bgs_solver.set_lhs_data (a);
      bgs_solver.print_blocks ();
      std::cout << "rank = " << rank << std::endl;
      toc ("set lhs data");
    }

  if (rank == 0)
    tic ();
  
  bgs_solver.factorize ();

  if (rank == 0)
    {
      std::cout << "rank = " << rank << std::endl;
      toc ("factorize");
    }
  
  if (rank == 0)
    {
      tic ();
      bgs_solver.set_rhs (rhs);
      std::cout << "rank = " << rank << std::endl;
      toc ("set rhs");
    }

  if (rank == 0)
    tic ();
  
  bgs_solver.solve ();

  if (rank == 0)
    {
      bgs_solver.print_resnorm ("first_solve_resnorm");
      std::cout << "rank = " << rank << std::endl;
      toc ("solve");
    }

  // if (rank == 0)
  //   for (int k = 0; k < rhs.size (); ++k)
  //     std::cout << "rhs[" << k << "] = " << rhs[k] << std::endl;
    
  // second solve
  if (rank == 0)
    {
      tic ();
      for (unsigned int k = 0; k < sp.rows (); ++k)
        sp[k][k] = 400.0;
      sp.aij_update (a, i, j);
      bgs_solver.set_lhs_data (a);
      std::cout << "rank = " << rank << std::endl;
      toc ("set lhs data");
    }

  //MPI_Barrier (MPI_COMM_WORLD);
  
  if (rank == 0)
    tic ();
  
  bgs_solver.factorize ();

  if (rank == 0)
    {
      std::cout << "rank = " << rank << std::endl;
      toc ("factorize");
    }

  if (rank == 0)
    {
      tic ();
      rhs.assign (rhs.size (), 2.0);
      bgs_solver.set_rhs (rhs);
      std::cout << "rank = " << rank << std::endl;
      toc ("set rhs");
    }  

  if (rank == 0)
    tic ();

  bgs_solver.solve ();

  if (rank == 0)
    {
      std::cout << "rank = " << rank << std::endl;
      bgs_solver.print_resnorm ("second_solve_resnorm");
      toc ("solve");
    }

  // if (rank == 0)
  //   for (int k = 0; k < rhs.size (); ++k)
  //     std::cout << "rhs[" << k << "] = " << rhs[k] << std::endl;

  if (rank == 0)
    tic ();
  
  bgs_solver.cleanup ();

  if (rank == 0)
    {
      std::cout << "rank = " << rank << std::endl;
      toc ("cleanup");
      print_timing_report ();
    }
  
  MPI_Finalize ();
  return (0);
}
