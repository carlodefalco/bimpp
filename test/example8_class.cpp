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
example8::set_rhs_values (const std::vector<double> & f_)
{
  f = f_;
};

void
example8::operator () (sparse_matrix& lhs,
                         std::vector<double>& rhs,
                         const std::vector<double>& guess)
{

  if (lhs.size() != n || rhs.size()!=n)
  {
     rhs.clear ();
     rhs.assign (n , 0.0);
     lhs.clear ();   // E' MEGLIO COSI' ?
     std::map<int, double> temp;
     temp.insert(std::pair<int,double>(0, 2*guess[0]) );
     lhs.push_back(temp);
     rhs[0]= - guess[0]*guess[0] + 1;      
     temp.clear();
     for (int i = 1; i<n-1 ; ++i)
     {  temp.insert (std::pair<int,double>(i-1, 1));
        temp.insert(std::pair<int,double>(i, -3*guess[i]*guess[i]));
        lhs.push_back(temp); 
        rhs[i]= - guess[i-1] + guess[i]*guess[i]*guess[i];
        temp.clear();
     }
     temp.insert(std::pair<int,double>(n-2,1));
     temp.insert (std::pair<int,double>(n-1,-1));
     lhs.push_back(temp); 
     rhs[n-1] = - guess[n-2] + guess[n-1];   
  }
  else
  {
     std::map<int, double> temp;
     temp.insert(std::pair<int,double>(0, 2*guess[0]) );
     lhs[0][0] = 2*guess[0];
     rhs[0]= - guess[0]*guess[0] + 1;      
     temp.clear();
     for (int i = 1; i<n-1 ; ++i)
     {  
        lhs[i][i-1] = 1;
        lhs[i][i] = -3*guess[i]*guess[i];
        rhs[i]= - guess[i-1] + guess[i]*guess[i]*guess[i];
     }
     lhs[n-1][n-2] = 1;
     lhs[n-1][n-1] = -1; 
     rhs[n-1] = - guess[n-2] + guess[n-1];

  }
}

void
example8::operator () (std::vector<double>& functional,
                        const std::vector<double>& guess)
{
  sparse_matrix M;
  operator () (M, functional, guess);
  for (unsigned int i = 0; i < functional.size (); ++i)
    functional[i] *= -1;
}

void
example8::get_exact_solution (std::vector<double> &exact_solution_)
{
  exact_solution_ = exact_solution;
}
