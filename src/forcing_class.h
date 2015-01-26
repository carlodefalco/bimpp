/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file forcing_class.h
  \brief interface for a forcing term.
*/

#ifndef HAVE_FORCING_CLASS_H
#define HAVE_FORCING_CLASS_H

#include "abstract_forcing_term.h"
#include "bim_sparse.h"
#include "operators.h"

/// Class that compute forcing term
/// \f$ \|F(x_k)-F(x_{k-1})-DF(x_{k-1})du_{k-1}\|/\|F(x_{k-1})\| \f$
class forcing_type1 : public abstract_forcing_term
{
private :
  double eta_max;

public :

  forcing_type1 (double eta_max_) :
    abstract_forcing_term ("Forcing Type 1"),
    eta_max (eta_max_) { };

  double
  operator () (const std::vector<double>& functional_old,
              const std::vector<double>& functional_new,
	      const std::vector<double>& df_gap,
              double eta_old);

};

/// Class that compute forcing term
/// \f$ \|F(x_k)\|-\|F(x_{k-1})+DF(x_{k-1})du_{k-1}\|/\|F(x_{k-1})\| \f$
class forcing_type2 : public abstract_forcing_term
{
private :
  double eta_max;

public :

  forcing_type2 (double eta_max_) :
    abstract_forcing_term ("Forcing Type 2"),
    eta_max (eta_max_) { };

  double
  operator () (const std::vector<double>& functional_old,
              const std::vector<double>& functional_new,
	      const std::vector<double>& df_gap,
              double eta_old);
};

/// Class that compute forcing term 
/// \f$ \gamma*(\|F(x_k)\|/\|F(x_{k-1})||)^{\alpha} \f$
class forcing_type3 : public abstract_forcing_term
{
private :

  double gamma;
  double alpha;
  double eta_max;

public :

  forcing_type3 (double gamma_, double alpha_, double eta_max_) :
    abstract_forcing_term ("Forcing Type 3"),
    gamma (gamma_),
    alpha (alpha_),
    eta_max (eta_max_) { };

  double
  operator () (const std::vector<double>& functional_old,
              const std::vector<double>& functional_new,
	      const std::vector<double>& df_gap,
              double eta_old);
};

/// Class that computes costant forcing term
class forcing_costant : public abstract_forcing_term
{
public :

  forcing_costant () :
    abstract_forcing_term ("Forcing Costant") { };

  double
  operator () (const std::vector<double>& functional_old,
              const std::vector<double>& functional_new,
	      const std::vector<double>& df_gap,
              double eta_old)
  { return eta_old; }
};

#endif
