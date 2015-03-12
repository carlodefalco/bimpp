/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file forcing_class.cpp
  \brief interfaces for a forcing term.
*/

#include "forcing_class.h"

double
forcing_type1::operator ()
(const std::vector<double>& f_old,
 const std::vector<double>& f_new,
 const std::vector<double>& df_gap,
 double eta_old)
{
  std::vector<double> temp (f_old.size (), 0.0);
  double eta_new = 0.0, eta_temp = 0.0;

  for (unsigned int i = 0; i < temp.size (); ++i)
    temp[i] = f_new[i] - f_old[i] - df_gap[i];

  for (unsigned int i = 0; i < temp.size (); ++i)
    {
      eta_new += temp[i] * temp[i];
      eta_temp += f_old[i] * f_old[i];
    }
  eta_new = sqrt (eta_new);
  eta_temp = sqrt (eta_temp);

  eta_new /= eta_temp;

  eta_temp = pow (eta_old, (1.0 + sqrt (5.0)) / 2.0);
  if (eta_temp > 0.1)
    return std::min (std::max (eta_new, eta_temp), eta_max);
  else return std::min (eta_new, eta_max);

}

double
forcing_type2::operator ()
(const std::vector<double>& f_old,
 const std::vector<double>& f_new,
 const std::vector<double>& df_gap,
 double eta_old)
{
  std::vector<double> temp (f_old.size (), 0.0);
  double eta_new = 0.0, eta_temp = 0.0;

  for (unsigned int i = 0; i < temp.size (); ++i)
    temp[i] = f_old[i] + df_gap[i];

  for (unsigned int i = 0; i < temp.size (); ++i)
    {
      eta_new += f_new[i] * f_new[i];
      eta_temp += temp[i] * temp[i];
    }
  eta_new = sqrt (eta_new);
  eta_temp = sqrt (eta_temp);

  eta_new -= eta_temp;

  eta_new = fabs (eta_new);

  eta_temp = 0.0;
  for (unsigned int i = 0; i < f_old.size (); ++i)
    eta_temp += f_old[i] * f_old[i];

  eta_temp = sqrt (eta_temp);

  eta_new /= eta_temp;

  eta_temp = pow (eta_old, (1.0 + sqrt (5.0)) / 2.0);
  if (eta_temp > 0.1)
    return std::min (std::max (eta_new, eta_temp), eta_max);
  else return std::min (eta_new, eta_max);
}

double
forcing_type3::operator ()
(const std::vector<double>& f_old,
 const std::vector<double>& f_new,
 const std::vector<double>& df_gap,
 double eta_old)
{
  std::vector<double> temp (f_old.size (), 0.0);
  double eta_new = 0.0, eta_temp = 0.0;

  for (unsigned int i = 0; i < temp.size (); ++i)
    {
      eta_new += f_new[i] * f_new[i];
      eta_temp += f_old[i] * f_old[i];
    }
  eta_new = sqrt (eta_new);
  eta_temp = sqrt (eta_temp);

  eta_new /= eta_temp;
  eta_new = gamma * pow (eta_new, alpha);
  eta_temp = gamma * pow (eta_old, alpha);

  if (eta_temp > 0.1)
    return std::min (std::max (eta_new, eta_temp), eta_max);
  else return std::min (eta_new, eta_max);

}
