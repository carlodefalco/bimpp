/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file plaplacian_class.cpp
  \brief interface for nonlinear problem
  \f$ -div (|\nabla u|^{p-2}\nabla u) = f \f$
*/

#include "plaplacian_class.h"

void
plaplacian::read_mesh (const std::string &mesh_name)
{
  msh.read (mesh_name);
  msh.precompute_properties ();
}

void
plaplacian::set_exact_solution
  (const std::vector<double> exact_solution_)
{
  exact_solution = exact_solution_;
};

void
plaplacian::set_rhs_values (const std::vector<double> & f_)
{
  f = f_;
};

void
plaplacian::set_boundary_conditions
  (std::vector<double> &boundary_values_,
   std::vector<int> &boundary_nodes_)
{
  boundary_values = boundary_values_;
  boundary_nodes = boundary_nodes_;
};

void
plaplacian::operator () (sparse_matrix& lhs,
                         std::vector<double>& rhs,
                         const std::vector<double>& guess)
{
  lhs.reset ();
  rhs.clear ();
  bim3a_structure (msh, lhs);

  std::vector<double> g;
  bim3a_pde_gradient (msh, guess, g);
  std::vector<double> modg (msh.nelements);

  for (int i = 0; i < msh.nelements; ++i)
    modg[i] = sqrt (g[3 * i + 0] * g[3 * i + 0] +
                    g[3 * i + 1] * g[3 * i + 1] +
                    g[3 * i + 2] * g[3 * i + 2]);

  std::vector<double> ecoeff (msh.nelements, 1.0);
  std::vector<double> isocoeff (msh.nelements, 1.0);
  std::vector<double> dcoeff (msh.nelements * 9, 1.0);

  for (int i = 0; i < msh.nelements; ++i)
    {
      dcoeff[0 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 0] * g[3 * i + 0];
      dcoeff[1 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 0] * g[3 * i + 1];
      dcoeff[2 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 0] * g[3 * i + 2];
      dcoeff[3 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 1] * g[3 * i + 0];
      dcoeff[4 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 1] * g[3 * i + 1];
      dcoeff[5 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 1] * g[3 * i + 2];
      dcoeff[6 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 2] * g[3 * i + 0];
      dcoeff[7 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 2] * g[3 * i + 1];
      dcoeff[8 + 9 * i] =
        pow (modg[i], p - 4.0) * (p - 2.0) *
        g[3 * i + 2] * g[3 * i + 2];

      isocoeff[i] = pow (modg[i], p - 2.0);
    }
  bim3a_laplacian (msh, isocoeff, lhs);
  rhs = lhs * guess;
  for (unsigned int i = 0; i < rhs.size (); ++i)
    rhs[i] *= -1;
  bim3a_rhs (msh, ecoeff, f, rhs);
  bim3a_laplacian_anisotropic (msh, dcoeff, lhs);
  bim3a_dirichlet_bc (lhs, rhs, boundary_nodes, boundary_values);
}

void
plaplacian::operator () (std::vector<double>& functional,
                        const std::vector<double>& guess)
{
  sparse_matrix M;
  operator () (M, functional, guess);
  for (unsigned int i = 0; i < functional.size (); ++i)
    functional[i] *= -1;
}

void
plaplacian::get_exact_solution (std::vector<double> &exact_solution_)
{
  exact_solution_ = exact_solution;
}
