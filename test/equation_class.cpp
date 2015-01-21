/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include "equation_class.h"

void
equation::operator () (sparse_matrix& lhs,
                       std::vector<double>& rhs,
                       const std::vector<double>& guess)
{
  lhs.resize (1);
  lhs[0][0] = 2.0;
  rhs.resize (1);
  rhs[0] = - (guess[0] * guess[0] - 2.0);
}

void
equation::operator () (std::vector<double>& functional,
                       const std::vector<double>& guess)
{
  functional.resize (1);
  functional[0] = guess[0] * guess[0] - 2.0;
}

void
equation::get_solution (std::vector<double>& sol)
{
  sol = exact_solution;
}
