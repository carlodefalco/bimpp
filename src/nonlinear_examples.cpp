/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <nonlinear_examples.h>

void
plaplacian::operator() (sparse_matrix& lhs,
                        std::vector<double>& rhs,
                        const std::vector<double>& u)
{
  lhs.reset ();
  rhs.clear ();
  bim3a_structure (msh, lhs);

  std::vector<double> g;
  bim3a_pde_gradient (msh, u, g);

  std::vector<double> modg (msh.nelements);
  for (int i = 0; i < this->msh.nelements; ++i)
    modg[i] = sqrt (g[3 * i + 0] * g[3 * i + 0] +
                    g[3 * i + 1] * g[3 * i + 1] +
                    g[3 * i + 2] * g[3 * i + 2]);

  std::vector<double> ecoeff (msh.nelements, 1.0);
  std::vector<double> isocoeff (msh.nelements, 1.0); //isotropic diffusion coefficients
  std::vector<double> dcoeff (msh.nelements * 9, 1.0);//anisotropic diffusion coefficients

  for (int i = 0; i < msh.nelements; ++i)
    {
      dcoeff[0 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 0] * g[3 * i + 0];
      dcoeff[1 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 0] * g[3 * i + 1];
      dcoeff[2 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 0] * g[3 * i + 2];
      dcoeff[3 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 1] * g[3 * i + 0];
      dcoeff[4 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 1] * g[3 * i + 1];
      dcoeff[5 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 1] * g[3 * i + 2];
      dcoeff[6 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 2] * g[3 * i + 0];
      dcoeff[7 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 2] * g[3 * i + 1];
      dcoeff[8 + 9 * i] = pow (modg[i], p - 4.0) * (p - 2.0) * g[3 * i + 2] * g[3 * i + 2];

      isocoeff[i] = pow (modg[i], p - 2.0);
    }
  bim3a_laplacian (msh, isocoeff, lhs);
  bim3a_matrix_vector_product (lhs, u, rhs);
  for (int i = 0; i < rhs.size (); ++i)
    rhs[i] *= -1;
  bim3a_rhs (msh, ecoeff, f, rhs);
  bim3a_laplacian_anisotropic (msh, dcoeff, lhs);
  bim3a_dirichletBC (lhs, rhs, bnodes, bsol);
}

void
plaplacian::operator() (std::vector<double>& F,const std::vector<double>& u)
{
  sparse_matrix M;
  this->operator()(M, F, u);
  for (int i = 0; i < F.size(); ++i)
    F[i] *= -1;
}

void
plaplacian::get_solution (std::vector<double>& sol)
{
  sol = exactsol;
}

void
equation::operator() (sparse_matrix& lhs,
                     std::vector<double>& rhs,
                     const std::vector<double>& u)
{
  lhs.resize (1);
  lhs[0][0] = 2.0;
  rhs.resize (1);
  rhs[0] = - (u[0] * u[0] - 2.0);
}

void
equation::operator() (std::vector<double>& f, const std::vector<double>& u)
{
  f.resize (1);
  f[0] = u[0] * u[0] - 2.0;
}

void
equation::get_solution (std::vector<double>& sol)
{
  sol = exactsol;
}
