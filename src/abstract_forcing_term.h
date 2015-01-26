/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file abstract_forcing_term.h
  \brief generic interface for a forcing term.
*/

#ifndef HAVE_ABSTRACT_FORCING_TERM_H
#define HAVE_ABSTRACT_FORCING_TERM_H 1

#include <string>
#include "abstract_nonlinear_problem.h"
#include "bim_sparse.h"
#include "operators.h"

class abstract_forcing_term
{
private :

  const std::string name;

protected :

  abstract_forcing_term (const char *name_) :
    name (name_) { };

public :

  /// Operator that returns the new forcing value
  /// for the nonlinear problem.
  /// Must be called on the master (rank == 0) node only.
  virtual double
    operator () (const std::vector<double>& functional_old,
                const std::vector<double>& functional_new,
                const std::vector<double>& df_gap,
                double old_forcing_value) = 0;

  /// Return the name of forcing term.
  const std::string&
  forcing_name () { return name; }

};

#endif
