/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file forcing_class.h
  \brief interfaces for a forcing term.
*/

#ifndef HAVE_FORCING_CLASS_H
#define HAVE_FORCING_CLASS_H

#include "abstract_forcing_term.h"
#include "bim_sparse.h"
#include "operators.h"

/// \brief Class that compute forcing term
/// \f$ ||F (x_k)-F (x_{k-1})-F' (x_{k-1})s_{k-1}||/||F (x_{k-1})||\f$.
/// \details safeguard:
/// \f$ \eta_k = max \{\eta_k, eta_{k-1}^{(1+\sqrt{5})/2)}\}\f$
/// if \f$  eta_{k-1}^{(1+\sqrt{5})/2)} > 0.1 \f$.
class forcing_type1 : public abstract_forcing_term
{
private :

  /// \brief Maximum value of forcing term.
  /// \details \f$ \eta_{max} \in (0, 1)\f$.
  double eta_max;

public :

  /// Default costructor.
  forcing_type1 (double eta_max_) :
    abstract_forcing_term ("Forcing Type 1"),
    eta_max (eta_max_) { };

  /// \brief Operator that returns the new forcing value
  /// for the nonlinear problem.
  double
  operator () (const std::vector<double>& functional_old,
               const std::vector<double>& functional_new,
               const std::vector<double>& df_gap,
               double eta_old);

};

/// \brief Class that compute forcing term
/// \f$ \big|||F (x_k)||-||F (x_{k-1})+F'(x_{k-1})s_{k-1}||\big|/||F (x_{k-1})||\f$.
/// \details safeguard:
/// \f$ \eta_k = max \{\eta_k, eta_{k-1}^{(1+\sqrt{5})/2)}\}\f$
/// if \f$  eta_{k-1}^{(1+\sqrt{5})/2)} > 0.1 \f$.
class forcing_type2 : public abstract_forcing_term
{
private :

  /// \brief Maximum value of forcing term.
  /// \details \f$ \eta_{max} \in (0, 1)\f$.
  double eta_max;

public :

  /// Default costructor
  forcing_type2 (double eta_max_) :
    abstract_forcing_term ("Forcing Type 2"),
    eta_max (eta_max_) { };

  /// \brief Operator that returns the new forcing value
  /// for the nonlinear problem.
  double
  operator () (const std::vector<double>& functional_old,
               const std::vector<double>& functional_new,
               const std::vector<double>& df_gap,
               double eta_old);
};

/// \brief Class that compute forcing term
/// \f$ \gamma (||F(x_k)||/||F(x_{k-1})||)^{\alpha}\f$.
/// \details safeguard:
/// \f$ \eta_k = max \{\eta_k, \gamma \eta_{k-1}^{\alpha}\}\f$
/// if \f$  \gamma \eta_{k-1}^{\alpha} > 0.1 \f$.
class forcing_type3 : public abstract_forcing_term
{
private :

  /// \brief Forcing parameter.
  /// \details \f$ \gamma \in [0, 1] \f$.
  double gamma;

  /// \brief Forcing parameter.
  /// \details \f$ \alpha \in (1, 2] \f$.
  double alpha;

  /// \brief Maximum vale of forcing term.
  /// \details \f$ \eta_{max} \in (0, 1)\f$.
  double eta_max;

public :

  /// Default costructor.
  forcing_type3 (double gamma_, double alpha_, double eta_max_) :
    abstract_forcing_term ("Forcing Type 3"),
    gamma (gamma_),
    alpha (alpha_),
    eta_max (eta_max_) { };

  /// \brief Operator that returns the new forcing value
  /// for the nonlinear problem.
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

  /// Default costructor
  forcing_costant () :
    abstract_forcing_term ("Forcing Costant") { };

  /// \brief Operator that returns the new forcing value
  /// for the nonlinear problem.
  double
  operator () (const std::vector<double>& functional_old,
               const std::vector<double>& functional_new,
               const std::vector<double>& df_gap,
               double eta_old)
  { return eta_old; }
};

#endif
