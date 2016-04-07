/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem:  
  \f[\frac{\partial u}{\partial t} -\Delta (u) = g\f]
  \f[u = 1-x^2-y^2-z^2 \:on \:boundary \f]
  \f[g = 6\f]

  Exact Solution:
  \f[u = 1-x^2-y^2-z^2\f]
*/

#include <stdio.h>
//#include <lis_config.h>
#include <lis.h>
#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <lis_class.h>
#include <fstream>
#include <stdlib.h>
#include <bim_config.h>

int main (int argc, char **argv)
{
  linear_solver *solver = new lis ();
  int rank, size;
  int base = solver->get_index_base ();

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  mesh msh;

  sparse_matrix       lhs, lhs_new;
  std::vector<double> rhs1, rhs2, rhs_new;
  std::vector<int>    ir, jc;
  std::vector<double> xa;
  std::vector<double> exactsolution_start, exactsolution;

  std::vector<int> bnodes;
  std::vector<double> vnodes_start, vnodes;

  std::vector<double> uold;
  double dt = 1;
  int T = 10;

  std::ofstream fout_sol ("Lis_Solution_TimeTest.txt");

  if (rank == 0)
    {
      std::cout << std::endl << std::endl
                << "*****\nLis test: Time Dependent Problem\n*****"
                << std::endl;

      std::cout << "read mesh" << std::endl;
      msh.read (data_dir + std::string ("mesh_in_cube.msh"));

      std::cout << "compute mesh props" << std::endl;
      msh.precompute_properties ();

      std::cout << "assemble stiffness matrix. nnodes = "
                << msh.nnodes << std::endl;

      bim3a_structure (msh, lhs);
       //isotropic diffusion coefficient
      std::vector<double> ecoeff (msh.nelements, 1.0);

      std::vector<double> v (msh.nnodes, 0.0);
      std::vector<double> ncoeff (msh.nnodes, 1 / dt);
      std::vector<double> nodecoeff1 (msh.nnodes, 1 / dt);
      std::vector<double> nodecoeff2 (msh.nnodes, 6.0);

      bim3a_advection_diffusion (msh, ecoeff, v, lhs);
      bim3a_reaction (msh, ecoeff, ncoeff, lhs);

      bim3a_rhs (msh, ecoeff, nodecoeff1, rhs1);
      bim3a_rhs (msh, ecoeff, nodecoeff2, rhs2);

      uold = std::vector<double> (msh.nnodes, 0.0);
      for (int i = 0; i < msh.nnodes; ++i)
        {
          uold[i] = 1.0 -
            msh.p (0, i) * msh.p (0, i) -
            msh.p (1, i) * msh.p (1, i) -
            msh.p (2, i) * msh.p (2, i);
        }

      std::vector<int> sidelist;
      sidelist.push_back (1);
      sidelist.push_back (2);
      sidelist.push_back (3);
      sidelist.push_back (4);
      sidelist.push_back (5);
      sidelist.push_back (6);

      msh.boundary_nodes (sidelist, bnodes);
      vnodes_start.resize (bnodes.size ());
      vnodes.resize (bnodes.size ());

      for (unsigned int i = 0; i < vnodes.size (); ++i)
        {
          vnodes_start[i] = 1.0 -
            msh.p (0, bnodes[i]) * msh.p (0, bnodes[i]) -
            msh.p (1, bnodes[i]) * msh.p (1, bnodes[i]) -
            msh.p (2, bnodes[i]) * msh.p (2, bnodes[i]);
        }

      exactsolution_start.resize (msh.nnodes);
      exactsolution.resize (msh.nnodes);

      for (unsigned int i = 0; i < exactsolution_start.size (); ++i)
        {
          exactsolution_start[i] = 1.0 -
            msh.p (0, i) * msh.p (0, i) -
            msh.p (1, i) * msh.p (1, i) -
            msh.p (2, i) * msh.p (2, i);
        }

      std::cout << std::endl << "Result of Time Dependent Test"
                << std::endl
                << "will be written in Lis_Solution_TimeTest.txt"
                << std::endl;
    }

  for (int t = 1; t <= T; ++t)
    {
      if (rank == 0)
        {
          std::cout << std::endl
                    << "Time Iteration: "<< t << std::endl;

          rhs_new.resize (rhs1.size ());
          lhs_new = lhs;

          for (unsigned int i = 0; i < rhs_new.size (); ++i)
            {
              rhs_new[i] = rhs2[i] + rhs1[i] * uold[i];
            }
          for (unsigned int i = 0; i < vnodes.size (); ++i)
            {
              vnodes[i] = vnodes_start[i];
            }

          bim3a_dirichlet_bc (lhs_new, rhs_new, bnodes, vnodes);

          for (unsigned int i = 0; i < exactsolution.size (); ++i)
            {
              exactsolution[i] = exactsolution_start[i];
            }
        }
      
      lhs_new.aij (xa, ir, jc, base);

     if (rank == 0)
        solver->set_lhs_structure (lhs_new.rows (), ir, jc);

     solver->analyze ();

      if (rank == 0)
        solver->set_lhs_data (xa);

     solver->factorize ();

      if (rank == 0)
        solver->set_rhs (rhs_new);
      if (solver->solver_type () == "iterative")
        {
          solver->set_max_iterations (1000);
          solver->set_tolerance (1e-12);
        }

     solver->solve ();

     if (rank == 0)
        {
          fout_sol << std::endl;
          fout_sol << "Time Iteration: "<< t << std::endl;

          double norm = 0;
          std::vector<double> delta (rhs_new.size ());

          for (unsigned int k = 0; k < rhs_new.size (); ++k)
            {
              fout_sol << rhs_new[k] << "  "
                       << exactsolution[k] << std::endl;
              delta[k] = exactsolution[k] - rhs_new[k];
            }
          fout_sol.close ();

          bim3a_norm (msh, delta, norm, Inf);
          std::cout << "Error: " << norm << std::endl;

          if (norm > 10e-10)
            {
              std::cerr << "The error is bigger than tolerance"
                        << std::endl;
              exit (-1);
            }
        }
    }
  fout_sol.close ();
  solver->cleanup ();
  return 0;
}
