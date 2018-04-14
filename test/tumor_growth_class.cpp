/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file tumor_growth_class.cpp
  
  \brief Interface for nonlinear problem
   \f[ 
      \begin{cases}
        \partial_t m - \mu \; div  (m \nabla p) = G(p) m \\
        \partial_t n - \nu \; div (n \nabla p) = 0 
      \end{cases}
   \f]
with \f$ p := K_{\gamma}(n+m)^{\gamma} \f$ , \f$ K_{\gamma} := \frac{\gamma + 1}{\gamma} \f$,

\f$ G(p) := \frac{200}{\pi} \arctan( 4 (p - P_M)) \f$ ,
 
\f$ m \f$ local density of dividing cells (tumor cells), 
\f$ n \f$ local density of non-dividing cells (not tumor cells) .

*/

#include "tumor_growth_class.h"

void
tumor_growth::set_initial_condition
  (const std::vector<double> & uold_)
{
  uold = uold_;
};

void
tumor_growth::set_rhs_values (const std::vector<double> & f_)
{
  f = f_;
};

void
tumor_growth::set_boundary_conditions
  (std::vector<double> &boundary_values_,
   std::vector<int> &boundary_nodes_)
{
  boundary_values = boundary_values_;
  boundary_nodes = boundary_nodes_;
};

void 
tumor_growth::set_matrices_structure ()
{

  mass.resize (n_nodes);
  Sm.resize (n_nodes);
  Sn.resize (n_nodes);
  mat_temp.resize (n_nodes);
  Amm.resize (n_nodes);
  Ann.resize (n_nodes);
  M.resize (2 * n_nodes);
  
  m.resize (n_nodes); 
  n.resize (n_nodes);
  mold.resize (n_nodes);
  nold.resize (n_nodes); 
  diffm.resize (n_nodes);
  diffn.resize (n_nodes);
  p.resize (n_nodes);
  G.resize (n_nodes);
  mdGdm.resize (n_nodes);
  mdGdn.resize (n_nodes);
  ecoeff.assign (n_elements, 1.0);
  ncoeff.assign (n_nodes, 1.0);
 
};

void
tumor_growth::operator () (sparse_matrix& lhs,
                         std::vector<double>& rhs,
                         const std::vector<double>& guess)
{
  if (rhs.size () != 2 * n_nodes)
    rhs.resize (2 * n_nodes);
  if (lhs.size () != 2 * n_nodes )
    lhs.resize (2 * n_nodes);
  
  lhs.reset ();
  
  Sm.reset ();
  Sn.reset ();
  mass.reset ();
  Amm.reset ();
  Ann.reset ();
  mat_temp.reset ();
 
  for (int i = 0; i < n_nodes; ++i)
    {
      m[i] = guess[i];
      mold[i] = uold[i];
      n[i] = guess[i+n_nodes];
      nold[i] = uold[i+n_nodes];
      
      diffm[i] = dt * mu * (g + 1) * m[i] * pow(n[i] + m[i], g-1); 
      diffn[i] = dt * nu * (g + 1) * n[i] * pow(n[i] + m[i], g-1);

      p[i] = (g + 1) / g * pow(n[i] + m[i], g);
      G[i] = ((200 / M_PI) * std::atan(4 * (PM - p[i])) >= 0 ? (200 / M_PI) * 
             std::atan(4 * (PM - p[i])) : 0);

      mdGdm[i] = (G[i] >= 0 ? m[i] * (200 / M_PI) * (- 4 * (g + 1) * pow(n[i] + m[i], g - 1))
                 / (1 + 16 * pow(PM - p[i], 2)) : 0);
      mdGdn[i] = (G[i] >= 0 ? m[i] * (200 / M_PI) * (- 4 * (g + 1) * pow(n[i] + m[i], g - 1)) 
                 / (1 + 16 * pow(PM - p[i], 2)) : 0);
     }

  bim2a_reaction (*tmsh, ecoeff, ncoeff, mass);
  
  bim2a_advection_eafe_diffusion (*tmsh, diffm, ncoeff, Sm); 
  bim2a_advection_eafe_diffusion (*tmsh, diffn, ncoeff, Sn);
  Sm.set_properties();
  Sn.set_properties();
  mass.set_properties();
  
  for (int i = 0; i < n_nodes; ++i)
    mat_temp[i][i] = - dt * mass[i][i] * G[i];
  mat_temp.set_properties();
 
  Amm = mass;
  Amm += Sm;
  Amm += mat_temp;
  Ann = mass;
  Ann += Sn;
  
  //rhs =  b - A * guess
 
  tempm = Amm * m;
  tempn = Sm * n;
  tempb = mass * mold;
  
  for (int i = 0; i < n_nodes; ++i)
    rhs[i] = - tempm[i] - tempn[i] + tempb[i] ;
    

  tempm = Sn * m;
  tempn = Ann * n;
  tempb = mass * nold;
  
  for (int i = 0; i <  n_nodes; ++i)
    rhs[i + n_nodes] = - tempm[i] - tempn[i] + tempb[i] ;
  
  //rhs1 =  -(Amm * m) - (Sm * n) + mass * mold ;
  //rhs2 =  -(Sn * m) - (Ann * n) + mass * nold ;
  
   
  mat_temp.reset ();
  for (int i = 0; i < n_nodes; ++i)
    mat_temp[i][i] = - dt * mass[i][i] * mdGdm[i];
  Amm += mat_temp;

  mat_temp.reset ();
  for (int i = 0; i < n_nodes; ++i)
    mat_temp[i][i] = - dt * mass[i][i] * mdGdn[i];
  Sm += mat_temp;
  
  for (int i = 0; i < n_nodes; ++i)
    for (auto j = Amm[i].begin (); j != Amm[i].end (); ++j)
      lhs[i][j->first] = j->second;

  for (int i = 0; i < n_nodes; ++i)
    for (auto j = Sm[i].begin (); j != Sm[i].end (); ++j)
      lhs[i][j->first + n_nodes] = j->second;
  
  for (int i = 0; i < n_nodes; ++i)
    for (auto j = Sn[i].begin (); j != Sn[i].end (); ++j)
      lhs[i + n_nodes][j->first] = j->second;
  
  for (int i = 0; i < n_nodes; ++i)
    for (auto j = Ann[i].begin (); j != Ann[i].end (); ++j)
      lhs[i + n_nodes][j->first + n_nodes] = j->second;

  
}


void
tumor_growth::operator () (std::vector<double>& functional,
                        const std::vector<double>& guess)
{
  //  sparse_matrix M;
  operator () (M, functional, guess);
  for (unsigned int i = 0; i < functional.size (); ++i)
    functional[i] *= -1;
}
