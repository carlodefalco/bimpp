/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <lis.h>
#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <lis_class.h>
#include <mpi.h>
#include <fstream>
#include <bim_config.h>

const int system_size = 10;
int shuffle (int x, int nnz)
{
  return (x < floor (nnz / 2) ?                                 
          x + floor (nnz / 2) + 1 :                             
          x - floor (nnz / 2));
}

int shuffle_not (int x, int nnz)
{ return (x); }

void
run_test_problem_rank0 (linear_solver *solver, int base, std::vector<double> &rhs,
                        int (*f) (int, int));

void
run_test_problem_rank1 (linear_solver *solver);

int main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  std::vector<double> lis_rhs, lis_rhs_shuffle;
  std::vector<double> mumps_rhs, mumps_rhs_shuffle;

  linear_solver *mumps_solver = new mumps ();
  int base = 1;

  if (rank == 0)
    run_test_problem_rank0 (mumps_solver, base,
                            mumps_rhs, shuffle_not);
  else
    run_test_problem_rank1 (mumps_solver);

  if (rank == 0)
    run_test_problem_rank0 (mumps_solver, base,
                            mumps_rhs_shuffle, shuffle);
  else
    run_test_problem_rank1 (mumps_solver);
  
 
  
  linear_solver *lis_solver = new lis ();
  lis_solver->set_tolerance (1e-12);
  lis_solver->set_preconditioner ("ilut");
  base = 0;

  if (rank == 0)
    run_test_problem_rank0 (lis_solver, base,
                            lis_rhs, shuffle_not);
  else
    run_test_problem_rank1 (lis_solver);

  if (rank == 0)
    run_test_problem_rank0 (lis_solver, base,
                            lis_rhs_shuffle, shuffle);
  else
    run_test_problem_rank1 (lis_solver);

  if (rank == 0)
    for (int ii = 0; ii < lis_rhs.size (); ++ii)
      std::cout << lis_rhs[ii]
                << lis_rhs_shuffle[ii]
                << "  " << mumps_rhs[ii]
                << "  " << mumps_rhs_shuffle[ii]
                << std::endl;
  
  mumps_solver->cleanup ();
  lis_solver->cleanup ();
  
  MPI_Finalize ();
  return (0);
}

void
run_test_problem_rank0 (linear_solver *solver, int base,
                        std::vector<double> &rhs, int (*f) (int, int))
{
  sparse_matrix       lhs;
  std::vector<int>    ir, jc, ir_tmp, jc_tmp;
  std::vector<double> xa, xa_tmp;

  lhs.resize (system_size);

  std::cout << "first solve using " << solver->solver_name () << std::endl;
  for (int ii = 0; ii < system_size; ++ii)
    {
      lhs[ii][ii] = 10;
      if (ii > 0)
        lhs[ii][ii-1] = -1;
      if (ii < system_size)
        lhs[ii][ii+1] = -1;
    }            

  // std::cout << lhs << std::endl;
  // std::cout << lhs.nnz << std::endl;
  
  std::cout << "\taij" << std::endl;
  lhs.aij (xa_tmp, ir_tmp, jc_tmp, base);
  xa = xa_tmp; ir = ir_tmp; jc = jc_tmp;
  for (int ii = 0; ii < xa_tmp.size (); ++ ii)
    {
      xa[ii] = xa_tmp[f (ii, xa_tmp.size ())];
      ir[ii] = ir_tmp[f (ii, ir_tmp.size ())];
      jc[ii] = jc_tmp[f (ii, jc_tmp.size ())];
    }
  // std::cout << xa.size () << std::endl;
  // std::cout << ir.size () << std::endl;
  // std::cout << jc.size () << std::endl;
  
  std::cout << "\tset_lhs_structure" << std::endl;
  solver->set_lhs_structure (lhs.rows (), ir, jc);

  std::cout << "\tanalyze" << std::endl;
  solver->analyze ();      

  std::cout << "\taij_update" << std::endl;
  lhs.aij_update (xa_tmp, ir_tmp, jc_tmp, base);
  for (int ii = 0; ii < xa_tmp.size (); ++ ii)
    xa[ii] = xa_tmp[f (ii, xa_tmp.size ())];

  for (int ii = 0; ii < xa.size (); ++ ii)
    std::cout << xa[ii] << std::endl;

  std::cout << "\tset_lhs_data" << std::endl;
  solver->set_lhs_data (xa);

  std::cout << "\tfactorize" << std::endl;
  solver->factorize ();

  std::cout << "\trhs.assign" << std::endl;
  rhs.assign (system_size, 1.0);

  std::cout << "\tset_rhs" << std::endl;
  solver->set_rhs (rhs);

  std::cout << "\tsolve" << std::endl;
  solver->solve ();

  std::cout << "second solve using " << solver->solver_name () << std::endl;
  for (int ii = 0; ii < system_size; ++ii)
    lhs[ii][ii] = 11;
  //  std::cout << lhs << std::endl;
  
  std::cout << "\taij_update" << std::endl;
  lhs.aij_update (xa_tmp, ir_tmp, jc_tmp, base);
  for (int ii = 0; ii < xa_tmp.size (); ++ ii)
    xa[ii] = xa_tmp[f (ii, xa_tmp.size ())];

  for (int ii = 0; ii < xa.size (); ++ ii)
    std::cout << xa[ii] << std::endl;
  
  std::cout << "\tset_lhs_data" << std::endl;
  solver->set_lhs_data (xa);
  
  std::cout << "\tfactorize" << std::endl;
  solver->factorize ();

  rhs.assign (system_size, 2.0);

  std::cout << "\tset_rhs" << std::endl;
  solver->set_rhs (rhs);

  std::cout << "\tsolve" << std::endl;
  solver->solve ();

};

void
run_test_problem_rank1 (linear_solver *solver)
{
  solver->analyze ();      


  solver->factorize ();
  solver->solve ();


  solver->factorize ();
  solver->solve ();

};
