/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file plaplacian_class.h
  \brief interface for nonlinear problem
  \f$ div (|\nabla u|^{p-2}\nabla u) = f \f$
*/

#ifndef HAVE_PLAPLACIAN_CLASS_H
#define HAVE_PLAPLACIAN_CLASS_H

#include "abstract_nonlinear_problem.h"
#include "operators.h"
#include "bim_sparse.h"

class plaplacian : public abstract_nonlinear_problem
{
private :

  double p;
  std::vector<double> exact_solution;
  std::vector<double> f;
  std::vector<double> boundary_values;
  std::vector<int> boundary_nodes;

public :

  plaplacian (double p_) :
    abstract_nonlinear_problem ("plaplacian"),
    p (p_) { };

  void
  read_mesh (const std::string &mesh_name);

  void
  set_exact_solution (const std::vector<double> exact_solution_);

  void
  set_rhs_values (const std::vector<double> & f_);

  void
  set_boundary_conditions
    (std::vector<double> &boundary_values_,
     std::vector<int> &boundary_nodes_);

  void
  operator () (sparse_matrix& lhs,
               std::vector<double>& rhs,
               const std::vector<double>& guess);
 
  void
  operator () (std::vector<double>& functional,
               const std::vector<double>& guess);

  void
  get_exact_solution (std::vector<double> &exact_solution_);
};

#endif
