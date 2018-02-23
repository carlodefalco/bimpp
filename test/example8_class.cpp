/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*!
  Problem:
  \f[ x_1^2 - 1 = 0 \f]
  \f[ x_1 - x_2^3 = 0\f]
  \f[ ... \f]
  \f[ x_{n-2} - x_{n-1}^3 = 0 \f]
  \f[ x_{n-1} - x_{n} = 0 \f]
 
  Exact Solution:
  \f[ x = (1, ... , 1)^T \f]
  
*/

#include "example8_class.h"

void
example8::set_exact_solution
  (const std::vector<double> exact_solution_)
{
  exact_solution = exact_solution_;
};


void
example8::operator () (sparse_matrix& lhs,
		       std::vector<double>& rhs,
		       const std::vector<double>& guess)
{

  rhs.resize (n);
  rhs.assign (n , 0.0);
  
  lhs.reset ();
  lhs.resize (n);

  
  
  lhs[0][0] = 2*guess[0];
  rhs[0]= - std::pow(guess[0],2) + 1;      

  
  for (int i = 1; i < n-1 ; ++i)
    {  
      lhs[i][i-1] = 1;
      lhs[i][i] = -3*std::pow(guess[i],2);
      rhs[i]= - guess[i-1] + std::pow(guess[i],3);
    }

  lhs[n-1][n-2] = 1;
  lhs[n-1][n-1] = -1; 
  rhs[n-1] = - guess[n-2] + guess[n-1];
  
}

void
example8::operator () (std::vector<double>& rhs,
		       const std::vector<double>& guess)
{


  rhs.resize (n);
  rhs.assign (n , 0.0);
  
  rhs[0]= -(- std::pow(guess[0],2) + 1);      

  for (int i = 1; i < n-1 ; ++i)
    rhs[i]= -(- guess[i-1] + std::pow (guess[i], 3));

  rhs[n-1] = -(- guess[n-2] + guess[n-1]);
}

void
example8::get_exact_solution (std::vector<double> &exact_solution_)
{
  exact_solution_ = exact_solution;
}
