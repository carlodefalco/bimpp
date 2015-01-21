/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file forcing_class.cpp
  \brief interface for a forcing term.
*/

#include "forcing_class.h"

double
forcing_type1::operator ()
(abstract_nonlinear_problem *problem,
 const std::vector<double>& old_guess,
 const std::vector<double>& gap_guess,
 double eta_old,
 norm_type norm_t)
{
  sparse_matrix df;
  std::vector<double> fnew, fold, unew, temp;
  double eta_new, eta_temp;
  unew.resize (gap_guess.size ());
  for (int i = 0; i < gap_guess.size (); ++i)
    unew[i] = old_guess[i] + gap_guess[i];

  (*problem) (fnew, unew);
  (*problem) (df, fold, old_guess);

  bim3a_matrix_vector_product (df, gap_guess, temp);

  for (unsigned int i = 0; i < temp.size (); ++i)
    temp[i] += fold[i] - fnew[i];

  bim3a_norm (problem->msh, fnew, eta_new, norm_t);
  bim3a_norm (problem->msh, fold, eta_temp, norm_t);

  eta_new /= eta_temp;

  eta_temp = pow (eta_old, (1.0 + sqrt (5.0)) / 2.0);
  if (eta_temp > 0.1)
    return std::min (std::max (eta_new, eta_temp), eta_max);
  else return std::min (eta_new, eta_max);
}

double
forcing_type2::operator ()
(abstract_nonlinear_problem *problem,
 const std::vector<double>& old_guess,
 const std::vector<double>& gap_guess,
 double eta_old,
 norm_type norm_t)
{
  sparse_matrix df;
  std::vector<double> fnew, fold, unew, temp;
  double eta_new, eta_temp;
  unew.resize (gap_guess.size ());
  for (int i = 0; i < gap_guess.size (); ++i)
    unew[i] = old_guess[i] + gap_guess[i];

  (*problem) (fnew, unew);
  (*problem) (df, fold, old_guess);

  bim3a_matrix_vector_product (df, gap_guess, temp);

  for (int i = 0;i < temp.size (); ++i)
    {
      temp[i]-=fold[i];
    }

  bim3a_norm (problem->msh, fnew, eta_new, norm_t);
  bim3a_norm (problem->msh, temp, eta_temp, norm_t);

  eta_new -= eta_temp;

  bim3a_norm (problem->msh, fold, eta_temp, norm_t);

  eta_new /= eta_temp;

  eta_temp = pow (eta_old, (1.0 + sqrt (5.0)) / 2.0);
  if (eta_temp > 0.1)
    return std::min (std::max (eta_new, eta_temp), eta_max);
  else return std::min (eta_new, eta_max);
}

double
forcing_type3::operator ()
(abstract_nonlinear_problem *problem,
 const std::vector<double>& old_guess,
 const std::vector<double>& gap_guess,
 double eta_old,
 norm_type norm_t)
{
  std::vector<double> fnew, fold, unew, temp;
  double eta_new, eta_temp;
  unew.resize (gap_guess.size ());
  for (int i = 0; i < gap_guess.size (); ++i)
    unew[i] = old_guess[i] + gap_guess[i];

  (*problem) (fnew, unew);
  (*problem) (fold, old_guess);
  bim3a_norm (problem->msh, fnew, eta_new, norm_t);
  bim3a_norm (problem->msh, fold, eta_temp, norm_t);

  eta_new /= eta_temp;
  eta_new = gamma * pow (eta_new, alpha);
  eta_temp = gamma * pow (eta_old, alpha);

  if (eta_temp > 0.1)
    return std::min (std::max (eta_new, eta_temp), eta_max);
  else return std::min (eta_new, eta_max);
}

