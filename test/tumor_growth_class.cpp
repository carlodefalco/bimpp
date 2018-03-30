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
tumor_growth::read_mesh (const std::string &mesh_name)
{
  const char* name = mesh_name.c_str();
  tmsh.load (name);
};

void
tumor_growth::set_exact_solution
  (const std::vector<double> & exact_solution_)
{
  exact_solution = exact_solution_;
};

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
tumor_growth::operator () (sparse_matrix& lhs,
                         std::vector<double>& rhs,
                         const std::vector<double>& guess)
{
  lhs.reset ();
  rhs.clear ();

  int n_nodes = tmsh.num_global_nodes (); 
  int n_elements = tmsh.num_local_quadrants ();
  
  sparse_matrix Amm, Amn, Anm, Ann, mass, Sm, Sn, matG, mGGm , mGGn;
  //bim2a_structure (msh, lhs_piece);

  std::vector<double> m (n_nodes), n (n_nodes), mold (n_nodes), nold (n_nodes);
  std::vector<double> diffm (n_nodes), diffn (n_nodes), p (n_nodes);
  std::vector<double> G (n_nodes), mdGdm (n_nodes), mdGdn (n_nodes);

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

  std::vector<double> ecoeff (n_elements, 1.0);
  std::vector<double> ncoeff (n_nodes, 1.0);

  mass.resize (n_nodes);
  bim2a_reaction (tmsh, ecoeff, ncoeff, mass);
  
  Sm.resize(n_nodes);
  Sn.resize(n_nodes);
  bim2a_advection_eafe_diffusion (tmsh, diffm, ncoeff, Sm); 
  bim2a_advection_eafe_diffusion (tmsh, diffn, ncoeff, Sn);
  
  matG.resize(n_nodes);
  mGGm.resize(n_nodes);
  mGGn.resize(n_nodes); 
  for (int i = 0; i < n_nodes; ++i)
    {
      matG[i][i] = - dt * mass[i][i] * G[i];
      mGGm[i][i] = - dt * mass[i][i] * mdGdm[i];
      mGGn[i][i] = - dt * mass[i][i] * mdGdn[i];
    }

  Amm = mass;
  Amm += Sm;
  Amm += matG;
  Amn = Sm;    // posso anche non creare Amn 
  Ann = mass;
  Ann += Sn;
  Anm = Sn;   // posso anche non creare Anm
  
  //rhs =  b - A * guess
  std::vector<double> rhs1(n_nodes), rhs2(n_nodes);
  std::vector<double> tempm , tempn, tempb;
  tempm = Amm * m;
  tempn = Amn * n;
  tempb = mass * mold;
  for (int i = 0; i < n_nodes; ++i)
  {
     rhs1[i] = - tempm[i] - tempn[i] + tempb[i] ;
  }  
  
  tempm = Anm * m;
  tempn = Ann * n;
  tempb = mass * nold;
  for (int i = 0; i < n_nodes; ++i)
  {
     rhs2[i] = - tempm[i] - tempn[i] + tempb[i] ;
  } 

  //rhs1 =  -(Amm * m) - (Amn * n) + mass * mold ;
  //rhs2 =  -(Anm * m) - (Ann * n) + mass * nold ;

  rhs.resize(2 * n_nodes);
  for (int i = 0 ; i < n_nodes ; ++i)
    {
      rhs[i] =  rhs1[i]; 
    }
  for (int i = 0 ; i < n_nodes ; ++i)
    {
      rhs[i + n_nodes] =  rhs2[i]; 
    }

 
  Amm += mGGm;
  Amn += mGGn;
  std::vector<double> xa;
  std::vector<int> jc, ir; 
  int indx = 0;
  
  lhs.resize(2*n_nodes);
  Amm.aij(xa, ir, jc, indx);
  for (int i = 0 ; i < xa.size(); ++i)
    {
      lhs[ir[i]][jc[i]] = xa[i];
    }
  
  Amn.aij(xa, ir, jc, indx );
  for (int i = 0 ; i < xa.size(); ++i)
    {
      lhs[ir[i]][jc[i] + n_nodes] = xa[i];
    }
 
  Anm.aij(xa, ir, jc, indx);
  for (int i = 0 ; i < xa.size(); ++i)
    {
      lhs[ir[i] + n_nodes][jc[i]] = xa[i];
    }
  
  Ann.aij(xa, ir, jc, indx);
  for (int i = 0 ; i < xa.size(); ++i)
    {
      lhs[ir[i] + n_nodes][jc[i] + n_nodes] = xa[i];
    }

}

void
tumor_growth::operator () (std::vector<double>& functional,
                        const std::vector<double>& guess)
{
  sparse_matrix M;
  operator () (M, functional, guess);
  for (unsigned int i = 0; i < functional.size (); ++i)
    functional[i] *= -1;
}

void
tumor_growth::get_exact_solution (std::vector<double> &exact_solution_)
{
  exact_solution_ = exact_solution;
}
