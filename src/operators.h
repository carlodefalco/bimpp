// Copyright (C) 2004-2011  Carlo de Falco
//
// This file is part of:
//     secs3d - A 3-D Drift--Diffusion Semiconductor Device Simulator 
//
//  secs3d is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  secs3d is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with secs3d; If not, see <http://www.gnu.org/licenses/>.
//
//  author: Carlo de Falco     <cdf _AT_ users.sourceforge.net>

#ifndef HAVE_OPERATORS_H
#define HAVE_OPERATORS_H 1
#include <mesh.h>
#include <bim_sparse.h>
#include <cmath>

//namespace bim
//{

/// Virtual class defining the interface 
/// of a FEM coefficient functor.
class 
coefficient_functor
{
public:
  virtual const double& operator() (const int) = 0;
};

/// Specialization to use a vector as a 
/// coefficient functor.
class 
coefficient_vector : public coefficient_functor
{
private:
  const std::vector<double> *v;
public:
  /// Constructor
  coefficient_vector (const std::vector<double> *_v) : v(_v) {};
  /// 
  const double& operator() (const int ii) { return (*v)[ii]; };
};

/// Allocate the structure for a FEM matrix over the mesh msh.
void bim3a_structure (const mesh& msh, sparse_matrix& A);

/// Assemble the rhs of a FEM problem.
void bim3a_rhs (const mesh& msh, const std::vector<double>& ecoeff, 
                const std::vector<double>& ncoeff, std::vector<double>& b);

/// Assemble the mass matrix of a FEM problem.
void bim3a_reaction (const mesh& msh, const std::vector<double>& ecoeff, 
                     const std::vector<double>& ncoeff, sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM diffusion problem.
void bim3a_laplacian (const mesh& msh, const std::vector<double>& epsilon, 
                      sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM advection-diffusion problem.
void bim3a_advection_diffusion (const mesh& msh, const std::vector<double>& epsilon, 
                                const std::vector<double>& phi, sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM diffusion problem (by the Orthogonal Subdomain Collocation method).
void bim3a_osc_laplacian (const mesh& msh, const std::vector<double>& epsilon, 
                          sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM advection-diffusion problem (by the Orthogonal Subdomain Collocation method).
void bim3a_osc_advection_diffusion (const mesh& msh, const std::vector<double>& epsilon,
                                    const std::vector<double>& phi, sparse_matrix& A);

/// Robustly compute B(x) = x / (exp(x) - 1).
void bimu_bernoulli (double x, double &bp, double &bn);
//}

#endif
