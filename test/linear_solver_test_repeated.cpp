/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <bim_sparse.h>
#include <cmath>
#include <mumps_class.h>
#include <lis_class.h>
#include <mpi.h>

const int system_size = 111371;
int shuffle (int x, int nnz)
{
  int half = nnz / 2;
  return (x >= half ?
          x - half :
          nnz - half + x);
}

int shuffle_not (int x, int nnz)
{ return (x); }

void
run_test_problem_rank0 (linear_solver *solver,
                        std::vector<double> &rhs,
                        int (*f) (int, int), bool b = true);

void
run_test_problem_rank1 (linear_solver *solver, bool b = true);

int main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);

  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  std::vector<double> lis_rhs (system_size, 0.0),
    lis_rhs_shuffle (system_size, 0.0);
  std::vector<double> mumps_rhs, mumps_rhs_shuffle;

  linear_solver *mumps_solver = new mumps ();

  if (rank == 0)
    run_test_problem_rank0 (mumps_solver,
                            mumps_rhs, shuffle_not);
  else
    run_test_problem_rank1 (mumps_solver);

  if (rank == 0)
    run_test_problem_rank0 (mumps_solver,
                            mumps_rhs_shuffle, shuffle);
  else
    run_test_problem_rank1 (mumps_solver);


  linear_solver *lis_solver = new lis ();
  
  lis_solver->set_tolerance (1e-14);
  lis_solver->set_max_iterations (500);
  lis_solver->set_preconditioner ("jacobi");

  if (rank == 0)
    run_test_problem_rank0 (lis_solver,
                            lis_rhs, shuffle_not);
  else
    run_test_problem_rank1 (lis_solver);

  if (rank == 0)
    run_test_problem_rank0 (lis_solver,
                            lis_rhs_shuffle, shuffle, false);
  else
    run_test_problem_rank1 (lis_solver, false);

  if (rank == 0 && false)
    for (unsigned int ii = 0; ii < lis_rhs.size (); ++ii)
      {
        std::cout << lis_rhs[ii]
                  << "  " << lis_rhs_shuffle[ii]
                  << "  " << mumps_rhs[ii]
                  << "  " << mumps_rhs_shuffle[ii]
                  << std::endl;
      }

  double lis_shuffle_to_shuffle_not = 0,
    lis_shuffle_to_mumps = 0,
    lis_shuffle_to_mumps_shuffle = 0;

  if (rank == 0)
    {
      for (unsigned int ii = 0; ii < lis_rhs.size (); ++ii)
        {
          double tmp = fabs (lis_rhs_shuffle[ii] - lis_rhs[ii]);
          lis_shuffle_to_shuffle_not =
            (lis_shuffle_to_shuffle_not < tmp) ?
            tmp : lis_shuffle_to_shuffle_not;

          tmp = fabs (lis_rhs_shuffle[ii] - mumps_rhs[ii]);
          lis_shuffle_to_mumps = (lis_shuffle_to_mumps < tmp) ?
            tmp : lis_shuffle_to_mumps;

          tmp = fabs (lis_rhs_shuffle[ii] - mumps_rhs_shuffle[ii]);
          lis_shuffle_to_mumps_shuffle =
            (lis_shuffle_to_mumps_shuffle < tmp) ?
            tmp : lis_shuffle_to_mumps_shuffle;
        }

      std::cout << "ls2l = " << lis_shuffle_to_shuffle_not
                << " ls2m = " << lis_shuffle_to_mumps
                << " ls2ms = " << lis_shuffle_to_mumps_shuffle
                << std::endl; 

      assert (lis_shuffle_to_shuffle_not < 1e-10);
      assert (lis_shuffle_to_mumps < 1e-10);
      assert (lis_shuffle_to_mumps_shuffle < 1e-10);      
    }
  MPI_Barrier (MPI_COMM_WORLD);
  
  mumps_solver->cleanup ();
  lis_solver->cleanup ();

  MPI_Finalize ();
  return (0);
}

void
run_test_problem_rank0 (linear_solver *solver,
                        std::vector<double> &rhs,
                        int (*f) (int, int), bool b)
{

  int base = solver->get_index_base ();

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
      if (ii < system_size - 1)
        lhs[ii][ii+1] = -1;
    }

  // std::cout << lhs << std::endl;
  // std::cout << lhs.nnz << std::endl;

  std::cout << "\taij" << std::endl;
  lhs.aij (xa_tmp, ir_tmp, jc_tmp, base);
  xa = xa_tmp; ir = ir_tmp; jc = jc_tmp;
  for (unsigned int ii = 0; ii < xa_tmp.size (); ++ ii)
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
  for (unsigned int ii = 0; ii < xa_tmp.size (); ++ ii)
    xa[ii] = xa_tmp[f (ii, xa_tmp.size ())];

  // for (int ii = 0; ii < xa.size (); ++ ii)
  //   std::cout << xa[ii] << std::endl;

  std::cout << "\tset_lhs_data" << std::endl;
  solver->set_lhs_data (xa);

  std::cout << "\tfactorize" << std::endl;
  solver->factorize ();

  std::cout << "\trhs.assign" << std::endl;
  rhs.assign (system_size, 3.0);

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
  for (unsigned int ii = 0; ii < xa_tmp.size (); ++ ii)
    xa[ii] = xa_tmp[f (ii, xa_tmp.size ())];

  // for (int ii = 0; ii < xa.size (); ++ ii)
  //   std::cout << xa[ii] << std::endl;

  std::cout << "\tset_lhs_data" << std::endl;
  solver->set_lhs_data (xa);

  std::cout << "\tfactorize" << std::endl;
  solver->factorize ();

  // std::cout.setf (std::ios::scientific, std::ios::floatfield);
  // std::cout.precision (17);
  
  for (int isolve = 0; isolve < 10; ++isolve)
    {
      std::cout << "\trhs.assign" << std::endl;
      rhs.assign (system_size, 2.0);
      
      std::cout << "\tsolve" << std::endl;
      solver->solve ();
      
      MPI_Barrier (MPI_COMM_WORLD);
    }
};

void
run_test_problem_rank1 (linear_solver *solver, bool b)
{

  // master node : solver->set_lhs_structure 
  solver->analyze ();
  // master node : solver->set_lhs_data
  
  solver->factorize ();
  // master node : solver->set_rhs
  
  solver->solve ();

  // master node : solver->set_lhs_data 
  solver->factorize ();

  for (int isolve = 0; isolve < 10; ++isolve)
    {
      solver->solve ();
      MPI_Barrier (MPI_COMM_WORLD);
    }

};
