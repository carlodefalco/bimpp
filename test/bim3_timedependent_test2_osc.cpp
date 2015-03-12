/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem:
  \f[ \frac{\partial u}{\partial t} - D\Delta (u) = g \f]
  \f[ u = (1-x^2-y^2-z^2)\cdot e^{-t} \:on \:boundary \f]
  \f[ g = (11+x^2+y^2+z^2)\cdot e^{-t} \f]
  \f[ D = diag (1.0, 2.0, 3.0) \f]

  Exact Solution:
  \f[ u = (1-x^2-y^2-z^2)\cdot e^{-t} \f]
*/

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <bim_config.h>

int main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  mesh msh;

  sparse_matrix       lhs;
  sparse_matrix       lhs_new;
  std::vector<double> rhs1;
  std::vector<double> rhs2;
  std::vector<double> rhs_new;
  std::vector<int>    ir, jc;
  std::vector<double> xa;

  std::vector<double> exactsolution_start, exactsolution;

  std::vector<int>    bnodes;
  std::vector<double> vnodes_start, vnodes;

  std::vector<double> uold;

  double dt = 0.2;
  int T = 5;
  char nomefile[] = "Solution_TimeTest2_osc.txt";
  std::ofstream fout (nomefile);

  for (int t = 1; t <= T; ++t)
    {
      if (rank == 0)
        {
          if (t == 1)
            {
              std::cout << std::endl << std::endl
                        << "*****" << std::endl
                        << "Time Dependent Test 2 (OSC)"
                        << "*****" << std::endl;

              std::cout << "read mesh" << std::endl;
              msh.read (data_dir + std::string ("mesh_in_cube.msh"));

              std::cout << "compute mesh props" << std::endl;
              msh.precompute_properties ();

              bim3a_structure (msh, lhs);
              std::vector<double> ecoeff (msh.nelements, 1.0);
              std::vector<double> dcoeff (msh.nelements * 3, 1.0);

              std::vector<double> v (msh.nnodes, 0.0);
              std::vector<double> ncoeff (msh.nnodes, 1 / dt);
              std::vector<double> nodecoeff1 (msh.nnodes, 1 / dt);
              std::vector<double> nodecoeff2 (msh.nnodes, 0.0);

              for (int i = 0; i < msh.nelements; ++i)
                {
                  dcoeff[0 + 3 * i] = 1.0;
                  dcoeff[1 + 3 * i] = 2.0;
                  dcoeff[2 + 3 * i] = 3.0;
                }

              for (int i = 0; i < msh.nnodes; ++i)
                nodecoeff2[i] = 12.0 - 1.0 +
                  msh.p (0, i) * msh.p (0, i) +
                  msh.p (1, i) * msh.p (1, i) +
                  msh.p (2, i) * msh.p (2, i);

              bim3a_osc_advection_diffusion_anisotropic
                (msh, dcoeff, v, lhs);
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

              for (unsigned int i = 0;
                   i < exactsolution_start.size (); ++i)
                {
                  exactsolution_start[i] = 1.0 -
                    msh.p (0, i) * msh.p (0, i) -
                    msh.p (1, i) * msh.p (1, i) -
                    msh.p (2, i) * msh.p (2, i);
                }
            }

          rhs_new.resize (rhs1.size ());
          lhs_new = lhs;

          for (unsigned int i = 0; i < rhs_new.size (); ++i)
            {
              rhs_new[i] = rhs2[i] * exp (- t * dt) +
                           rhs1[i] * uold[i];
            }
          for (unsigned int i = 0;i < vnodes.size (); ++i)
            {
              vnodes[i] = vnodes_start[i] * exp (- t * dt);
            }

          bim3a_dirichlet_bc (lhs_new, rhs_new, bnodes, vnodes);

          if (t == 1)
            lhs_new.aij (xa, ir, jc, 1);
          else
            lhs_new.aij_update (xa, ir, jc, 1);

          for (unsigned int i = 0; i < exactsolution.size (); ++i)
            {
              exactsolution[i] =
                exactsolution_start[i] * exp (- t * dt);
            }
        }

      mumps mumps_solver;

      if (rank == 0)
        mumps_solver.set_lhs_structure (lhs_new.rows (), ir, jc);

      mumps_solver.analyze ();

      if (rank == 0)
        mumps_solver.set_lhs_data (xa);

      mumps_solver.factorize ();

      if (rank == 0)
        mumps_solver.set_rhs (rhs_new);

      mumps_solver.solve ();

      if (rank == 0)
        {
          if (t == 1)
            std::cout << std::endl << "Result of Time Dependent Test"
                      << std::endl << "will be written in "
                      << nomefile << std::endl;

          std::cout << std::endl << "Iteration " << t << std::endl;
          fout << std::endl << "Iteration " << t << std::endl;

          double norm = 0;
          std::vector<double> delta (rhs_new.size ());

          for (unsigned int k = 0; k < rhs_new.size (); ++k)
            {
              uold[k] = rhs_new[k];
              fout << rhs_new[k] << "  "
                   << exactsolution[k] << std::endl;
              delta[k] = exactsolution[k] - rhs_new[k];
            }

          fout.close ();

          bim3a_norm (msh, delta, norm, Inf);
          std::cout << "Error: " << norm << std::endl;

          if (norm > 10e-3)
            {
              std::cerr << "The error is bigger than tolerance"
                        << std::endl;
              exit (-1);
            }
        }

      mumps_solver.cleanup ();

    }
  fout.close ();
  MPI_Finalize ();

  return (0);
}
