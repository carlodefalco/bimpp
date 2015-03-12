/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem: 
  \f[ -div (D (\nabla (u)-\nabla (v)u)) + u = g \f]
  \f[ u = 1-2x^2-2y^2-z^2 \:on \:boundary \f]
  \f[ g = 7-2x^2-2y^2-z^2-2x-2y-2z \f]
  \f[ D = diag (0.5,0.5,1) \f]
  \f[ v = x + y + z \f]

  Exact Solution:
  \f[ u=1-2x^2-2y^2-z^2 \f]
*/

#include <bim_sparse.h>
#include <mesh.h>
#include <operators.h>
#include <mumps_class.h>
#include <mpi.h>
#include <fstream>
#include <bim_config.h>

int main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  mesh msh;

  sparse_matrix       lhs;
  std::vector<double> rhs;
  std::vector<int>    ir, jc;
  std::vector<double> xa;
  std::vector<double> exactsolution;
  if (rank == 0)
    {
      std::cout << "\n\n*****\nStationary Test 3 (OSC)\n*****\n";

      std::cout << "read mesh" << std::endl;
      msh.read (data_dir + std::string ("mesh_in_cube.msh"));

      std::cout << "compute mesh props" << std::endl;
      msh.precompute_properties ();

      bim3a_structure (msh, lhs);
      //anisotropic diffusion coefficient
      std::vector<double> dcoeff (msh.nelements * 3, 1.0);

      std::vector<double> ecoeff (msh.nelements, 1.0);
      std::vector<double> v (msh.nnodes, 0.0);
      std::vector<double> ncoeff (msh.nnodes, 0.0);
      std::vector<double> nodecoeff (msh.nnodes,1.0);

      for (int k = 0; k < msh.nelements; ++k)
        {
          dcoeff[0 + 3 * k] = 0.5;
          dcoeff[1 + 3 * k] = 0.5;
          dcoeff[2 + 3 * k] = 1;
        }

      for (int k = 0; k < msh.nnodes; ++k)
        {
          v[k] = msh.p (0, k) + msh.p (1, k) + msh.p (2, k);
          ncoeff[k] = 7.0 -
            2 * msh.p (0, k) * msh.p (0, k) -
            2 * msh.p (1, k) * msh.p (1, k) -
            msh.p (2, k) * msh.p (2, k) -
            2 * msh.p (0, k) - 2 * msh.p (1, k) - 2 * msh.p (2, k);
        }

      bim3a_osc_advection_diffusion_anisotropic (msh, dcoeff, v, lhs);

      bim3a_reaction (msh, ecoeff, nodecoeff, lhs);
      bim3a_rhs (msh, ecoeff, ncoeff, rhs);

      std::vector<int> sidelist;
      sidelist.push_back (1);
      sidelist.push_back (2);
      sidelist.push_back (3);
      sidelist.push_back (4);
      sidelist.push_back (5);
      sidelist.push_back (6);

      std::vector<int> bnodes;
      std::vector<double> vnodes;

      msh.boundary_nodes (sidelist, bnodes);
      vnodes.resize (bnodes.size ());

      for (unsigned int i = 0; i < vnodes.size (); ++i)
        {
          vnodes[i] = 1.0 -
            2 * msh.p (0, bnodes[i]) * msh.p (0, bnodes[i]) -
            2 * msh.p (1, bnodes[i]) * msh.p (1, bnodes[i]) -
            msh.p (2, bnodes[i]) * msh.p (2, bnodes[i]);
        }

      bim3a_dirichlet_bc (lhs, rhs, bnodes, vnodes);

      lhs.aij (xa, ir, jc, 1);

      exactsolution.resize (msh.nnodes);

      for (unsigned int i = 0; i < exactsolution.size (); ++i)
        {
          exactsolution[i] = 1.0 -
            2 * msh.p (0, i) * msh.p (0, i) -
            2 * msh.p (1, i) * msh.p (1, i) -
            msh.p (2, i) * msh.p (2, i);
        }
    }

  mumps mumps_solver;

  if (rank == 0)
    mumps_solver.set_lhs_structure (lhs.rows (), ir, jc);

  mumps_solver.analyze ();

  if (rank == 0)
    mumps_solver.set_lhs_data (xa);

  mumps_solver.factorize ();

  if (rank == 0)
    mumps_solver.set_rhs (rhs);

  mumps_solver.solve ();

  if (rank == 0)
    {
      std::cout << "\nResult of Stationary Test"
                << std::endl
                << "will be written in solution3_OSC.txt"
                << std::endl;
      std::ofstream fout ("solution3_OSC.txt");
      fout << std::endl;

      double norm = 0;
      std::vector<double> delta (rhs.size ());

      for (unsigned int k = 0; k < rhs.size (); ++k)
        {
          fout << rhs[k] << "  " << exactsolution[k] << std::endl;
          delta[k] = exactsolution[k] - rhs[k];
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

  MPI_Finalize ();
  return (0);
}
