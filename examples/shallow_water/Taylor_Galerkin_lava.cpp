#include "Taylor_Galerkin_lava.h"
#include <algorithm>
#include <cassert>

TG2_scheme::TG2_scheme(Q1& sol,
                       Q1& sold,
                       Q1& soldd, 
                       Q1& incr,
                       std::vector<std::array<double,4>>& incr_anti_diff,
                       Q1& P_plus,
                       Q1& P_minus,
                       Q0& sol_onehalf,
                       Q1& mass,
                       const ordering& oh,
                       const ordering& oUx,
                       const ordering& oUy,
                       const ordering& oTh,
                       const Q1& Z,
                       Q0& Z_onehalf,
                       const double& DELTAT,
                       const double& h_min,
                       const bool& is_non_reflBC,
                       const bool& is_isothermal,
                       const double& grav,
                       const double& nu_ref,
                       const double& T_ref,
                       const double& b_exp_coeff,
                       const double& density)
: sol(sol), sold(sold), soldd(soldd), incr(incr), incr_anti_diff(incr_anti_diff), P_plus(P_plus), P_minus(P_minus), sol_onehalf(sol_onehalf), mass(mass), 
  ordh(oh), ordUx(oUx), ordUy(oUy), ordTh(oTh), Z(Z), Z_onehalf(Z_onehalf), DELTAT(DELTAT), epsilon(h_min), is_non_reflBC(is_non_reflBC), grav(grav), nu_ref(nu_ref), T_ref(T_ref),
  density(density), is_isothermal(is_isothermal), b_exp_coeff(b_exp_coeff)
{ }
 
 

void
TG2_scheme::compute_dt (tmesh::quadrant_iterator quadrant)
{
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii); 
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hdof[ii]  = sol [ordh  (quadrant->gt (ii) )];
      Uxdof[ii] = sol [ordUx (quadrant->gt (ii) )];
      Uydof[ii] = sol [ordUy (quadrant->gt (ii) )]; 
    }
    else
    {
      hdof[ii]  = .5 * (sol [ordh  (quadrant->gparent (0, ii) )] +
                        sol [ordh  (quadrant->gparent (1, ii) )]);
      Uxdof[ii] = .5 * (sol [ordUx (quadrant->gparent (0, ii) )] +
                        sol [ordUx (quadrant->gparent (1, ii) )]);
      Uydof[ii] = .5 * (sol [ordUy (quadrant->gparent (0, ii) )] +
                        sol [ordUy (quadrant->gparent (1, ii) )]);
    }
    const double vel_abs = hdof[ii]>epsilon ? std::sqrt(Uxdof [ii]*Uxdof [ii] + Uydof [ii]*Uydof [ii])/hdof[ii] : 0.; 
    //Fr = 10.;//.*(hdof[ii]>epsilon); //(vel_abs>1.e14 ? 1. : 0.);
    //Fr = vel_abs>std::sqrt(grav*hdof[ii]) ? 0. : 0.;
    Fr = std::max(hdof[ii]>epsilon ? vel_abs/std::sqrt(grav*hdof[ii]) : 0., Fr);

    //if (Fr!=0) std::cout << Fr << std::endl;
  }
  
  for (int ii = 0; ii < 4; ++ii){
    const auto& hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    
    const auto vel_rusanov_cell_x = hpoint>epsilon ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
    const auto vel_rusanov_cell_y = hpoint>epsilon ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;
    
    const auto dtoptx = hpoint>epsilon ? Dx/vel_rusanov_cell_x : DELTAT;
    const auto dtopty = hpoint>epsilon ? Dy/vel_rusanov_cell_y : DELTAT;
    const auto dtopt = dtoptx > dtopty ? dtopty : dtoptx;

    if (dt > dtopt) set_dt (dtopt);
    
  }
}

void
TG2_scheme::compute_dt_adaptive (tmesh::quadrant_iterator quadrant)
{
  
  double Nu_hmean_cell = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    double hdof, hdof_old, hdof_oldd, dh_t, h1, h2, h3, a_coeff, b_coeff;
    if (! quadrant->is_hanging (ii) )
    {
      hdof      = sol   [ordh  (quadrant->gt (ii) )];
      hdof_old  = sold  [ordh  (quadrant->gt (ii) )];
      hdof_oldd = soldd [ordh  (quadrant->gt (ii) )];
    }
    else
    {
      hdof      = .5 * (sol   [ordh  (quadrant->gparent (0, ii) )] +
                        sol   [ordh  (quadrant->gparent (1, ii) )]);
      hdof_old  = .5 * (sold  [ordh  (quadrant->gparent (0, ii) )] +
                        sold  [ordh  (quadrant->gparent (1, ii) )]);
      hdof_oldd = .5 * (soldd [ordh  (quadrant->gparent (0, ii) )] +
                        soldd [ordh  (quadrant->gparent (1, ii) )]);
    }
    dh_t = (hdof - hdof_old)/(time - timed);
    
    h1   = hdof_oldd/((timedd - timed )*(timedd - time ));
    h2   = hdof_old /((timed  - timedd)*(timed  - time ));
    h3   = hdof     /((time   - timedd)*(time   - timed));
    
    a_coeff = h1+h2+h3;
    b_coeff = - (h1*(time+timed) + h2*(time+timedd) + h3*(timed+timedd));
    
    Nu_hmean_cell += (1./3.*a_coeff*a_coeff*(time*time+time*timed+timed*timed) + a_coeff*(b_coeff-dh_t)*(time+timed) + (b_coeff-dh_t)*(b_coeff-dh_t));
    
  }
  Nu_hmean_cell /= 4.;
  nu_htot += Nu_hmean_cell*(time-timed)*(time-timed); // the dimension is length^2, this is eta^2
  
}

void
TG2_scheme::first_step (tmesh::quadrant_iterator quadrant)
{
  
  const auto & index_quadrant_global = quadrant->get_global_quad_idx (); 
  
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii);
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;


  
  double h_cell_average = 0., Ux_cell_average = 0., Uy_cell_average = 0., Th_cell_average = 0.;// source_Ux_cell_average = 0., source_Uy_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    double hdof_c, Uxdof_c, Uydof_c, Thdof_c;

    if (! quadrant->is_hanging (ii) )
    {
      hdof_c    = sol [ordh    (quadrant->gt (ii) )];
      Uxdof_c   = sol [ordUx   (quadrant->gt (ii) )];
      Uydof_c   = sol [ordUy   (quadrant->gt (ii) )];
      Thdof_c   = sol [ordTh   (quadrant->gt (ii) )];

      Z_node[ii]  = Z [quadrant->gt (ii)];
    }
    else
    {
      hdof_c    = .5 * (sol [ordh    (quadrant->gparent (0, ii) )] +
                        sol [ordh    (quadrant->gparent (1, ii) )]);
      Uxdof_c   = .5 * (sol [ordUx   (quadrant->gparent (0, ii) )] +
                        sol [ordUx   (quadrant->gparent (1, ii) )]);
      Uydof_c   = .5 * (sol [ordUy   (quadrant->gparent (0, ii) )] +
                        sol [ordUy   (quadrant->gparent (1, ii) )]);
      Thdof_c   = .5 * (sol [ordTh   (quadrant->gparent (0, ii) )] +
                        sol [ordTh   (quadrant->gparent (1, ii) )]);

      Z_node[ii]  = .5 * (Z [quadrant->gparent(0,ii)] +
                          Z [quadrant->gparent(1,ii)]);
    }

    h_cell_average  += hdof_c;
    Ux_cell_average += Uxdof_c;
    Uy_cell_average += Uydof_c;
    Th_cell_average += Thdof_c;

    
    fluxx_h_node[ii]  = h_flux_formula_x   (hdof_c, Uxdof_c, Uydof_c);
    fluxy_h_node[ii]  = h_flux_formula_y   (hdof_c, Uxdof_c, Uydof_c);
    fluxx_Ux_node[ii] = Ux_flux_formula_x  (hdof_c, Uxdof_c, Uydof_c);
    fluxy_Ux_node[ii] = Ux_flux_formula_y  (hdof_c, Uxdof_c, Uydof_c);
    fluxx_Uy_node[ii] = Uy_flux_formula_x  (hdof_c, Uxdof_c, Uydof_c);
    fluxy_Uy_node[ii] = Uy_flux_formula_y  (hdof_c, Uxdof_c, Uydof_c);
    fluxx_Th_node[ii] = Th_flux_formula_x  (hdof_c, Uxdof_c, Uydof_c, Thdof_c);
    fluxy_Th_node[ii] = Th_flux_formula_y  (hdof_c, Uxdof_c, Uydof_c, Thdof_c);

  }
  h_cell_average    /= 4.;
  Ux_cell_average   /= 4.;
  Uy_cell_average   /= 4.;
  Th_cell_average   /= 4.;
  
  //source_Ux_cell_average /= 4.;
  //source_Uy_cell_average /= 4.;

  
  const auto div_Fh_x = .5*((fluxx_h_node[1]-fluxx_h_node[0]) + (fluxx_h_node[3]-fluxx_h_node[2]));
  const auto div_Fh_y = .5*((fluxy_h_node[2]-fluxy_h_node[0]) + (fluxy_h_node[3]-fluxy_h_node[1]));
  const auto div_Fh_cell = Dy*div_Fh_x + Dx*div_Fh_y;
  
  const auto div_FUx_x = .5*((fluxx_Ux_node[1]-fluxx_Ux_node[0]) + (fluxx_Ux_node[3]-fluxx_Ux_node[2]));
  const auto div_FUx_y = .5*((fluxy_Ux_node[2]-fluxy_Ux_node[0]) + (fluxy_Ux_node[3]-fluxy_Ux_node[1]));
  const auto div_FUx_cell = Dy*div_FUx_x + Dx*div_FUx_y;
  
  const auto div_FUy_x = .5*((fluxx_Uy_node[1]-fluxx_Uy_node[0]) + (fluxx_Uy_node[3]-fluxx_Uy_node[2]));
  const auto div_FUy_y = .5*((fluxy_Uy_node[2]-fluxy_Uy_node[0]) + (fluxy_Uy_node[3]-fluxy_Uy_node[1]));
  const auto div_FUy_cell = Dy*div_FUy_x + Dx*div_FUy_y;

  const auto div_FTh_x = .5*((fluxx_Th_node[1]-fluxx_Th_node[0]) + (fluxx_Th_node[3]-fluxx_Th_node[2]));
  const auto div_FTh_y = .5*((fluxy_Th_node[2]-fluxy_Th_node[0]) + (fluxy_Th_node[3]-fluxy_Th_node[1]));
  const auto div_FTh_cell = Dy*div_FTh_x + Dx*div_FTh_y;

  const auto slope_x_c = ((Z_node[1] - Z_node[0]) + (Z_node[3] - Z_node[2]))/Dx/2.;
  const auto slope_y_c = ((Z_node[2] - Z_node[0]) + (Z_node[3] - Z_node[1]))/Dy/2.;


  Z_onehalf[index_quadrant_global] = (Z_node[0]+Z_node[1]+Z_node[2]+Z_node[3])*.25;
  
  sol_onehalf[ordh    (index_quadrant_global)] = h_cell_average  - dt*.5 *  div_Fh_cell /area;
  sol_onehalf[ordUx   (index_quadrant_global)] = Ux_cell_average - dt*.5 * (div_FUx_cell/area - src_slope_formula (h_cell_average, slope_x_c));
  sol_onehalf[ordUy   (index_quadrant_global)] = Uy_cell_average - dt*.5 * (div_FUy_cell/area - src_slope_formula (h_cell_average, slope_y_c));
  sol_onehalf[ordTh   (index_quadrant_global)] = Th_cell_average - dt*.5 *  div_FTh_cell/area;
  
  // add source term for the momentum
  const auto & h_onehalf_updated = sol_onehalf[ordh    (index_quadrant_global)];

  sol_onehalf[ordUx   (index_quadrant_global)] /= 1. - dt*.5*Ux_src_formula(h_onehalf_updated, 1., 0., Th_cell_average);
  sol_onehalf[ordUy   (index_quadrant_global)] /= 1. - dt*.5*Uy_src_formula(h_onehalf_updated, 0., 1., Th_cell_average);
  
} 


void
TG2_scheme::compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant)
{

  // look at tmesh.h
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 
  const auto & index_quadrant_global = quadrant->get_global_quad_idx (); 

  std::array<int,4> bimpp_to_rev_ord = {0, 1, 3, 2};
  
  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;

  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){
      hdof  [ii]  = sol [ordh    (quadrant->gt (ii))];
      Uxdof [ii]  = sol [ordUx   (quadrant->gt (ii))];
      Uydof [ii]  = sol [ordUy   (quadrant->gt (ii))];
      Thdof [ii]  = sol [ordTh   (quadrant->gt (ii))];

      Z_node[ii]  = Z [quadrant->gt (ii)];
      
      isdof_or_hanging[ii] = 1.;
    } else {
      hdof  [ii]  = .5 * (sol [ordh  (quadrant->gparent(0,ii))] +
                          sol [ordh  (quadrant->gparent(1,ii))]);
      Uxdof [ii]  = .5 * (sol [ordUx (quadrant->gparent(0,ii))] +
                          sol [ordUx (quadrant->gparent(1,ii))]);
      Uydof [ii]  = .5 * (sol [ordUy (quadrant->gparent(0,ii))] +
                          sol [ordUy (quadrant->gparent(1,ii))]);
      Thdof [ii]  = .5 * (sol [ordTh (quadrant->gparent(0,ii))] +
                          sol [ordTh (quadrant->gparent(1,ii))]);

      Z_node[ii]  = .5 * (Z [quadrant->gparent(0,ii)] +
                          Z [quadrant->gparent(1,ii)]);
      
      isdof_or_hanging[ii] = .5;
    }
  }

  // weights coefficients for the flux term
  der_coeffs_x = {-Dy/2.*isdof_or_hanging[0], +Dy/2.*isdof_or_hanging[1],
    -Dy/2.*isdof_or_hanging[2], +Dy/2.*isdof_or_hanging[3]};
   
  der_coeffs_y = {-Dx/2.*isdof_or_hanging[0], -Dx/2.*isdof_or_hanging[1],
    +Dx/2.*isdof_or_hanging[2], +Dx/2.*isdof_or_hanging[3]};

  double vel_rusanov_cell_x = 0., vel_rusanov_cell_y = 0.;
  for (int ii = 0; ii < 4; ++ii){
    const auto& hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    vel_rusanov_cell_x += hpoint>epsilon ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
    vel_rusanov_cell_y += hpoint>epsilon ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.; 

  const std::array<double,4> eta_vec = {hdof[0]+Z_node[0]*g_coeff, hdof[1]+Z_node[1]*g_coeff, hdof[2]+Z_node[2]*g_coeff, hdof[3]+Z_node[3]*g_coeff};

  grad_cell_Z    = {.5 * ( (Z_node [3] - Z_node [2]) + (Z_node [1] - Z_node [0]) ), .5 * ( (Z_node [2] - Z_node [0]) + (Z_node [3] - Z_node [1]) )};
  grad_cell_eta  = {.5 * ( (eta_vec[3] - eta_vec[2]) + (eta_vec[1] - eta_vec[0]) ), .5 * ( (eta_vec[2] - eta_vec[0]) + (eta_vec[3] - eta_vec[1]) )};
  grad_cell_h    = {.5 * ( (hdof   [3] - hdof   [2]) + (hdof   [1] - hdof   [0]) ), .5 * ( (hdof   [2] - hdof   [0]) + (hdof   [3] - hdof   [1]) )};
  grad_cell_Ux   = {.5 * ( (Uxdof  [3] - Uxdof  [2]) + (Uxdof  [1] - Uxdof  [0]) ), .5 * ( (Uxdof  [2] - Uxdof  [0]) + (Uxdof  [3] - Uxdof  [1]) )};
  grad_cell_Uy   = {.5 * ( (Uydof  [3] - Uydof  [2]) + (Uydof  [1] - Uydof  [0]) ), .5 * ( (Uydof  [2] - Uydof  [0]) + (Uydof  [3] - Uydof  [1]) )};
  grad_cell_Th   = {.5 * ( (Thdof  [3] - Thdof  [2]) + (Thdof  [1] - Thdof  [0]) ), .5 * ( (Thdof  [2] - Thdof  [0]) + (Thdof  [3] - Thdof  [1]) )};
 
  //std::cout << grad_cell_eta[0] << " " << grad_cell_eta[1] << " " << grad_cell_eta[2] << " " << grad_cell_eta[3] << std::endl;

  const double & h_cell    = sol_onehalf[ordh    (index_quadrant_global)];
  const double & Ux_cell   = sol_onehalf[ordUx   (index_quadrant_global)];
  const double & Uy_cell   = sol_onehalf[ordUy   (index_quadrant_global)];
  const double & Th_cell   = sol_onehalf[ordTh   (index_quadrant_global)];

  const double & Z_cell    = Z_onehalf[index_quadrant_global];

  const auto slope_x_cell = 0;//grad_cell_Z[0];
  const auto slope_y_cell = 0;//grad_cell_Z[1];

  //std::cout << slope_x_cell/Dx << " " << slope_y_cell/Dy << std::endl;

  const auto contr_slope_x = .5*src_slope_formula (h_cell, slope_x_cell)*Dy;
  const auto contr_slope_y = .5*src_slope_formula (h_cell, slope_y_cell)*Dx;


  std::array<double,4> contr_x = {0., 0., 0., 0.}, contr_y = {0., 0., 0., 0.};


  // boundary conditions, just for the transport term!
  bool is_boundary_edge = true;
  for (int iEdge = 0; iEdge < 4; ++iEdge){

    is_boundary_edge = true;

    const auto i_1 = bimpp_to_rev_ord[iEdge];
    const auto i_2 = bimpp_to_rev_ord[(iEdge+1)%4];

    const auto edge_length = std::sqrt(std::pow((xn[i_1]-xn[i_2]),2.) + std::pow((yn[i_1]-yn[i_2]),2.));
    const std::array<double,2> outward_normal_edge = {(-yn[i_1]+yn[i_2])/edge_length, ( xn[i_1]-xn[i_2])/edge_length};  

    for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
         quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
    {
      std::array<double,4> Xn, Yn;

      for (int ii = 0; ii < 4; ++ii) {
        Xn[ii] = quadrant_nei->p(0, ii);
        Yn[ii] = quadrant_nei->p(1, ii);

        if (! quadrant_nei->is_hanging (ii)){
          Z_node_nei[ii]  = Z [quadrant_nei->gt (ii)];
        } else {
          Z_node_nei[ii]  = .5 * (Z [quadrant_nei->gparent(0,ii)] +
                                  Z [quadrant_nei->gparent(1,ii)]);
        }
  
      }

      const auto & index_quadrant_nei_global = quadrant_nei->get_global_quad_idx (); 

      for (int jEdge = 0; jEdge < 4; ++jEdge) { // cycle neigh edges 

        const auto j_1 = bimpp_to_rev_ord[jEdge];
        const auto j_2 = bimpp_to_rev_ord[(jEdge+1)%4];

        const auto edge_length_nei = std::sqrt(std::pow((Xn[j_1]-Xn[j_2]),2.) + std::pow((Yn[j_1]-Yn[j_2]),2.));
        const std::array<double,2> outward_normal_edge_nei = {(-Yn[j_1]+Yn[j_2])/edge_length_nei, ( Xn[j_1]-Xn[j_2])/edge_length_nei};  
        const bool check_orthogonality = std::inner_product(outward_normal_edge_nei.begin(), outward_normal_edge_nei.end(), outward_normal_edge.begin(), 0.) == -1;

        bool is_owned_quadrant = index_quadrant_nei_global>=Z_onehalf.get_range_start () && index_quadrant_nei_global<Z_onehalf.get_range_end ();

        if ( (((xn[i_1] == Xn[j_1] && yn[i_1] == Yn[j_1]) || 
               (xn[i_2] == Xn[j_1] && yn[i_2] == Yn[j_1]))||
              ((xn[i_1] == Xn[j_2] && yn[i_1] == Yn[j_2]) ||
               (xn[i_2] == Xn[j_2] && yn[i_2] == Yn[j_2]))) && check_orthogonality && index_quadrant_global!=index_quadrant_nei_global )
        {
          is_boundary_edge = false;

          // scrivere qui la somma dei contributi per i termini non-cons.!
          double Z_cell_nei = is_owned_quadrant ? Z_onehalf[index_quadrant_nei_global] : (Z_node_nei[0]+Z_node_nei[1]+Z_node_nei[2]+Z_node_nei[3])*.25;
          
          // .5 salta fuori dall'integrazione per trapezi tra 0 e 1 in coordinata \xi (è il valore in LHS da metter qui sotto!)
          contr_x[i_1] += .5*signum(outward_normal_edge[0])*grav*h_cell*(Z_cell_nei - Z_cell)*isdof_or_hanging[i_1];
          contr_x[i_2] += .5*signum(outward_normal_edge[0])*grav*h_cell*(Z_cell_nei - Z_cell)*isdof_or_hanging[i_2];

          contr_y[i_1] += .5*signum(outward_normal_edge[1])*grav*h_cell*(Z_cell_nei - Z_cell)*isdof_or_hanging[i_1];
          contr_y[i_2] += .5*signum(outward_normal_edge[1])*grav*h_cell*(Z_cell_nei - Z_cell)*isdof_or_hanging[i_2];

          //break; // this just goes outside the jEdge cycle 
        }
      }

    }

    if (is_boundary_edge) // set boundary conditions
    { 
      /*
      contr_x[i_1] += grav*h_cell*slope_x_cell*isdof_or_hanging[i_1]*Dx;
      contr_x[i_2] += grav*h_cell*slope_x_cell*isdof_or_hanging[i_2]*Dx;

      contr_y[i_1] += grav*h_cell*slope_y_cell*isdof_or_hanging[i_1]*Dy;
      contr_y[i_2] += grav*h_cell*slope_y_cell*isdof_or_hanging[i_2]*Dy;
      */
      

      auto h_cell_nei  = h_cell;
      auto Ux_cell_nei = Ux_cell;
      auto Uy_cell_nei = Uy_cell;
      auto Th_cell_nei = Th_cell;

      Ux_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell + outward_normal_edge[1]*Uy_cell)*outward_normal_edge[0];
      Uy_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell + outward_normal_edge[1]*Uy_cell)*outward_normal_edge[1];

      const auto speed     = h_cell    >epsilon ? std::abs((Ux_cell    /h_cell    )*outward_normal_edge[0] + (Uy_cell    /h_cell    )*outward_normal_edge[1]) + std::sqrt(grav*h_cell    ) : 0.;
      const auto speed_nei = h_cell_nei>epsilon ? std::abs((Ux_cell_nei/h_cell_nei)*outward_normal_edge[0] + (Uy_cell_nei/h_cell_nei)*outward_normal_edge[1]) + std::sqrt(grav*h_cell_nei) : 0.;

      const auto smax = std::max(speed, speed_nei); 

      const auto flux_int_h  = .5*((h_flux_formula_x (h_cell, Ux_cell, Uy_cell)+h_flux_formula_x (h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (h_flux_formula_y (h_cell, Ux_cell, Uy_cell)+h_flux_formula_y (h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(h_cell_nei -h_cell );
      const auto flux_int_Ux = .5*((Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(Ux_cell_nei-Ux_cell);
      const auto flux_int_Uy = .5*((Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uy_cell_nei-Uy_cell);
      const auto flux_int_Th = .5*((Th_flux_formula_x(h_cell, Ux_cell, Uy_cell, Th_cell)+Th_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei, Th_cell_nei))*outward_normal_edge[0] + (Th_flux_formula_y(h_cell, Ux_cell, Uy_cell, Th_cell)+Th_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei, Th_cell_nei))*outward_normal_edge[1]) - .5*smax*(Th_cell_nei-Th_cell);

      //std::cout << flux_int_Ux << " " << Ux_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei) << " " << Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell) << std::endl;

      auto outward_normal_edge_abs = outward_normal_edge;
      outward_normal_edge_abs[0] = std::abs(outward_normal_edge[0]);
      outward_normal_edge_abs[1] = std::abs(outward_normal_edge[1]);

/*
      const auto diff_term_h_x  = grad_cell_h [0] * vel_rusanov_cell_y;
      const auto diff_term_h_y  = grad_cell_h [1] * vel_rusanov_cell_x;

      const auto diff_term_Ux_x = (is_non_reflBC ? grad_cell_Ux[0] : -grad_cell_Ux[0]) * vel_rusanov_cell_y;
      const auto diff_term_Ux_y = (is_non_reflBC ? grad_cell_Ux[1] : -grad_cell_Ux[1]) * vel_rusanov_cell_x;

      const auto diff_term_Uy_x = (is_non_reflBC ? grad_cell_Uy[0] : -grad_cell_Uy[0]) * vel_rusanov_cell_y;
      const auto diff_term_Uy_y = (is_non_reflBC ? grad_cell_Uy[1] : -grad_cell_Uy[1]) * vel_rusanov_cell_x;


      const auto flux_int_h  = .5*((h_flux_formula_x (h_cell, Ux_cell, Uy_cell)+h_flux_formula_x (h_cell_nei, Ux_cell_nei, Uy_cell_nei)-2.*diff_term_h_x               )*outward_normal_edge[0] + (h_flux_formula_y (h_cell, Ux_cell, Uy_cell)+h_flux_formula_y (h_cell_nei, Ux_cell_nei, Uy_cell_nei)-2.*diff_term_h_y               )*outward_normal_edge[1]);
      const auto flux_int_Ux = .5*((Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei)-2.*diff_term_Ux_x*is_non_reflBC)*outward_normal_edge[0] + (Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei)-2.*diff_term_Ux_y*is_non_reflBC)*outward_normal_edge[1]);
      const auto flux_int_Uy = .5*((Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei)-2.*diff_term_Uy_x*is_non_reflBC)*outward_normal_edge[0] + (Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei)-2.*diff_term_Uy_y*is_non_reflBC)*outward_normal_edge[1]);
*/

      
      // .5 is the base function evaluated in the middle, mid-point intergration
      if (! quadrant->is_hanging (i_1))
      {
        incr[ordh   (quadrant->gt (i_1))] += -edge_length*flux_int_h *.5;
        incr[ordUx  (quadrant->gt (i_1))] += -edge_length*flux_int_Ux*.5 + contr_slope_x*outward_normal_edge_abs[0];
        incr[ordUy  (quadrant->gt (i_1))] += -edge_length*flux_int_Uy*.5 + contr_slope_y*outward_normal_edge_abs[1];
        incr[ordTh  (quadrant->gt (i_1))] += -edge_length*flux_int_Th*.5;
      }
      else
      {
        incr [ordh  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_h *.5*.5;
        incr [ordh  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_h *.5*.5;
      
        incr [ordUx (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Ux*.5*.5 + contr_slope_x*.5*outward_normal_edge_abs[0];
        incr [ordUx (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Ux*.5*.5 + contr_slope_x*.5*outward_normal_edge_abs[0];
      
        incr [ordUy (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uy*.5*.5 + contr_slope_y*.5*outward_normal_edge_abs[1];
        incr [ordUy (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uy*.5*.5 + contr_slope_y*.5*outward_normal_edge_abs[1];

        incr [ordTh (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Th *.5*.5;
        incr [ordTh (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Th *.5*.5;
      }

      if (! quadrant->is_hanging (i_2))
      {
        incr[ordh   (quadrant->gt (i_2))] += -edge_length*flux_int_h *.5;
        incr[ordUx  (quadrant->gt (i_2))] += -edge_length*flux_int_Ux*.5 + contr_slope_x*outward_normal_edge_abs[0];
        incr[ordUy  (quadrant->gt (i_2))] += -edge_length*flux_int_Uy*.5 + contr_slope_y*outward_normal_edge_abs[1];
        incr[ordTh  (quadrant->gt (i_2))] += -edge_length*flux_int_Th *.5;
      }
      else
      {
        // il secondo .5 è per hanging nodes
        incr [ordh  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_h *.5*.5;
        incr [ordh  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_h *.5*.5;
      
        incr [ordUx (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Ux*.5*.5 + contr_slope_x*.5*outward_normal_edge_abs[0];
        incr [ordUx (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Ux*.5*.5 + contr_slope_x*.5*outward_normal_edge_abs[0];
      
        incr [ordUy (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uy*.5*.5 + contr_slope_y*.5*outward_normal_edge_abs[1];
        incr [ordUy (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uy*.5*.5 + contr_slope_y*.5*outward_normal_edge_abs[1];

        incr [ordTh (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Th *.5*.5;
        incr [ordTh (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Th *.5*.5;
      }

      /*
      auto h_cell_nei  = h_cell;
      auto Ux_cell_nei = Ux_cell;
      auto Uy_cell_nei = Uy_cell;

      Ux_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell + outward_normal_edge[1]*Uy_cell)*outward_normal_edge[0];
      Uy_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell + outward_normal_edge[1]*Uy_cell)*outward_normal_edge[1];

      const auto diff_term_h_x  = grad_cell_h [0] * vel_rusanov_cell_y;
      const auto diff_term_h_y  = grad_cell_h [1] * vel_rusanov_cell_x;

      const auto diff_term_Ux_x = (is_non_reflBC ? grad_cell_Ux[0] : -grad_cell_Ux[0]) * vel_rusanov_cell_y;
      const auto diff_term_Ux_y = (is_non_reflBC ? grad_cell_Ux[1] : -grad_cell_Ux[1]) * vel_rusanov_cell_x;

      const auto diff_term_Uy_x = (is_non_reflBC ? grad_cell_Uy[0] : -grad_cell_Uy[0]) * vel_rusanov_cell_y;
      const auto diff_term_Uy_y = (is_non_reflBC ? grad_cell_Uy[1] : -grad_cell_Uy[1]) * vel_rusanov_cell_x;


      const auto F_star_h_x  = h_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei) - diff_term_h_x;
      const auto F_star_h_y  = h_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei) - diff_term_h_y;

      const auto F_star_Ux_x = Ux_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei) - diff_term_Ux_x;
      const auto F_star_Ux_y = Ux_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei) - diff_term_Ux_y;

      const auto F_star_Uy_x = Uy_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei) - diff_term_Uy_x;
      const auto F_star_Uy_y = Uy_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei) - diff_term_Uy_y;

      const auto h_1    = der_coeffs_x[i_1]*F_star_h_x +der_coeffs_y[i_1]*F_star_h_y;
      const auto Ux_1   = der_coeffs_x[i_1]*F_star_Ux_x+der_coeffs_y[i_1]*F_star_Ux_y;
      const auto Uy_1   = der_coeffs_x[i_1]*F_star_Uy_x+der_coeffs_y[i_1]*F_star_Uy_y;

      const auto h_2    = der_coeffs_x[i_2]*F_star_h_x +der_coeffs_y[i_2]*F_star_h_y;
      const auto Ux_2   = der_coeffs_x[i_2]*F_star_Ux_x+der_coeffs_y[i_2]*F_star_Ux_y;
      const auto Uy_2   = der_coeffs_x[i_2]*F_star_Uy_x+der_coeffs_y[i_2]*F_star_Uy_y;

      if (! quadrant->is_hanging (i_1))
      {
        incr[ordh   (quadrant->gt (i_1))] += h_1;
        incr[ordUx  (quadrant->gt (i_1))] += Ux_1;
        incr[ordUy  (quadrant->gt (i_1))] += Uy_1;
      }
      else
      {
        incr [ordh  (quadrant->gparent(0,i_1))] += h_1;
        incr [ordh  (quadrant->gparent(1,i_1))] += h_1;
      
        incr [ordUx (quadrant->gparent(0,i_1))] += Ux_1;
        incr [ordUx (quadrant->gparent(1,i_1))] += Ux_1;
      
        incr [ordUy (quadrant->gparent(0,i_1))] += Uy_1;
        incr [ordUy (quadrant->gparent(1,i_1))] += Uy_1;
      }

      if (! quadrant->is_hanging (i_2))
      {
        incr[ordh   (quadrant->gt (i_2))] += h_2;
        incr[ordUx  (quadrant->gt (i_2))] += Ux_2;
        incr[ordUy  (quadrant->gt (i_2))] += Uy_2;
      }
      else
      {
        // il secondo .5 è per hanging nodes
        incr [ordh  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_h *.5*.5;
        incr [ordh  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_h *.5*.5;
      
        incr [ordUx (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Ux*.5*.5;
        incr [ordUx (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Ux*.5*.5;
      
        incr [ordUy (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uy*.5*.5;
        incr [ordUy (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uy*.5*.5;
      }*/

    }

  }
  

  const auto diff_term_h_x  = h_cell>epsilon ? grad_cell_eta[0]*vel_rusanov_cell_y*.5 : grad_cell_h[0]*vel_rusanov_cell_y*.5;
  const auto diff_term_h_y  = h_cell>epsilon ? grad_cell_eta[1]*vel_rusanov_cell_x*.5 : grad_cell_h[1]*vel_rusanov_cell_x*.5;

  //const auto diff_term_h_x  = grad_cell_eta[0]*vel_rusanov_cell_y*.5;
  //const auto diff_term_h_y  = grad_cell_eta[1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Ux_x = grad_cell_Ux [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Ux_y = grad_cell_Ux [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uy_x = grad_cell_Uy [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uy_y = grad_cell_Uy [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Th_x = grad_cell_Th [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Th_y = grad_cell_Th [1]*vel_rusanov_cell_x*.5;




  const auto F_star_h_x  = h_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_h_x;
  const auto F_star_h_y  = h_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_h_y;

  const auto F_star_Ux_x = Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_x;
  const auto F_star_Ux_y = Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_y;

  const auto F_star_Uy_x = Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_x;
  const auto F_star_Uy_y = Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_y;

  const auto F_star_Th_x = Th_flux_formula_x(h_cell, Ux_cell, Uy_cell, Th_cell) - diff_term_Th_x;
  const auto F_star_Th_y = Th_flux_formula_y(h_cell, Ux_cell, Uy_cell, Th_cell) - diff_term_Th_y;


  for (int ii = 0; ii < 4; ++ii){

    const auto h_    = der_coeffs_x[ii]*F_star_h_x +der_coeffs_y[ii]*F_star_h_y;
    const auto Ux_   = der_coeffs_x[ii]*F_star_Ux_x+der_coeffs_y[ii]*F_star_Ux_y - .5*Dy*contr_x[ii];// + .25*area*src_slope_formula(h_cell, slope_x_cell);
    const auto Uy_   = der_coeffs_x[ii]*F_star_Uy_x+der_coeffs_y[ii]*F_star_Uy_y - .5*Dx*contr_y[ii];// + .25*area*src_slope_formula(h_cell, slope_y_cell);
    const auto Th_   = der_coeffs_x[ii]*F_star_Th_x+der_coeffs_y[ii]*F_star_Th_y;

    const auto h_al  = der_coeffs_x[ii]*diff_term_h_x  + der_coeffs_y[ii]*diff_term_h_y; 
    const auto Ux_al = der_coeffs_x[ii]*diff_term_Ux_x + der_coeffs_y[ii]*diff_term_Ux_y;
    const auto Uy_al = der_coeffs_x[ii]*diff_term_Uy_x + der_coeffs_y[ii]*diff_term_Uy_y;
    const auto Th_al = der_coeffs_x[ii]*diff_term_Th_x + der_coeffs_y[ii]*diff_term_Th_y;

    incr_anti_diff[ordh (index_quadrant)][ii] = h_al;
    incr_anti_diff[ordUx(index_quadrant)][ii] = Ux_al;
    incr_anti_diff[ordUy(index_quadrant)][ii] = Uy_al;
    incr_anti_diff[ordTh(index_quadrant)][ii] = Th_al;


    if (! quadrant->is_hanging (ii)){

      incr [ordh  (quadrant->gt (ii))] += h_;
      incr [ordUx (quadrant->gt (ii))] += Ux_;
      incr [ordUy (quadrant->gt (ii))] += Uy_;
      incr [ordTh (quadrant->gt (ii))] += Th_;

      P_plus [ordh  (quadrant->gt (ii))] += std::max(0., h_al );
      P_plus [ordUx (quadrant->gt (ii))] += std::max(0., Ux_al);
      P_plus [ordUy (quadrant->gt (ii))] += std::max(0., Uy_al);
      P_plus [ordTh (quadrant->gt (ii))] += std::max(0., Th_al);

      P_minus [ordh  (quadrant->gt (ii))] += std::min(0., h_al );
      P_minus [ordUx (quadrant->gt (ii))] += std::min(0., Ux_al);
      P_minus [ordUy (quadrant->gt (ii))] += std::min(0., Uy_al);
      P_minus [ordTh (quadrant->gt (ii))] += std::min(0., Th_al);
      
      
    } else {


      incr [ordh  (quadrant->gparent(0,ii))] += h_;
      incr [ordh  (quadrant->gparent(1,ii))] += h_;
      
      incr [ordUx (quadrant->gparent(0,ii))] += Ux_;
      incr [ordUx (quadrant->gparent(1,ii))] += Ux_;
      
      incr [ordUy (quadrant->gparent(0,ii))] += Uy_;
      incr [ordUy (quadrant->gparent(1,ii))] += Uy_;

      incr [ordTh (quadrant->gparent(0,ii))] += Th_;
      incr [ordTh (quadrant->gparent(1,ii))] += Th_;



      P_plus [ordh  (quadrant->gparent(0,ii))] += std::max(0., h_al);
      P_plus [ordh  (quadrant->gparent(1,ii))] += std::max(0., h_al);
      
      P_plus [ordUx (quadrant->gparent(0,ii))] += std::max(0., Ux_al);
      P_plus [ordUx (quadrant->gparent(1,ii))] += std::max(0., Ux_al);
      
      P_plus [ordUy (quadrant->gparent(0,ii))] += std::max(0., Uy_al);
      P_plus [ordUy (quadrant->gparent(1,ii))] += std::max(0., Uy_al);

      P_plus [ordTh (quadrant->gparent(0,ii))] += std::max(0., Th_al);
      P_plus [ordTh (quadrant->gparent(1,ii))] += std::max(0., Th_al);



      P_minus [ordh  (quadrant->gparent(0,ii))] += std::min(0., h_al);
      P_minus [ordh  (quadrant->gparent(1,ii))] += std::min(0., h_al);
      
      P_minus [ordUx (quadrant->gparent(0,ii))] += std::min(0., Ux_al);
      P_minus [ordUx (quadrant->gparent(1,ii))] += std::min(0., Ux_al);
      
      P_minus [ordUy (quadrant->gparent(0,ii))] += std::min(0., Uy_al);
      P_minus [ordUy (quadrant->gparent(1,ii))] += std::min(0., Uy_al);

      P_minus [ordTh (quadrant->gparent(0,ii))] += std::min(0., Th_al);
      P_minus [ordTh (quadrant->gparent(1,ii))] += std::min(0., Th_al);
      
    }
  }



}



void
TG2_scheme::second_step (tmesh::quadrant_iterator quadrant)
{

  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 

  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }

  
  for (int ii = 0; ii < 4; ++ii){

    double hdof_c, Uxdof_c, Uydof_c, Thdof_c, P_plus_h_c, P_minus_h_c, P_plus_Ux_c, P_minus_Ux_c, P_plus_Uy_c, P_minus_Uy_c, P_plus_Th_c, P_minus_Th_c;

    if (! quadrant->is_hanging (ii)){
      hdof_c      = sol [ordh    (quadrant->gt (ii))];
      Uxdof_c     = sol [ordUx   (quadrant->gt (ii))];
      Uydof_c     = sol [ordUy   (quadrant->gt (ii))];
      Thdof_c     = sol [ordTh   (quadrant->gt (ii))];

      Z_node[ii]    = Z [quadrant->gt (ii)];

      P_plus_h_c   = P_plus [ordh   (quadrant->gt (ii))];
      P_minus_h_c  = P_minus[ordh   (quadrant->gt (ii))];

      P_plus_Ux_c   = P_plus [ordUx   (quadrant->gt (ii))];
      P_minus_Ux_c  = P_minus[ordUx   (quadrant->gt (ii))];

      P_plus_Uy_c   = P_plus [ordUy   (quadrant->gt (ii))];
      P_minus_Uy_c  = P_minus[ordUy   (quadrant->gt (ii))];

      P_plus_Th_c   = P_plus [ordTh   (quadrant->gt (ii))];
      P_minus_Th_c  = P_minus[ordTh   (quadrant->gt (ii))];
      

    } else {
      hdof_c   = .5 * (sol [ordh  (quadrant->gparent(0,ii))] +
                       sol [ordh  (quadrant->gparent(1,ii))]);
      Uxdof_c  = .5 * (sol [ordUx (quadrant->gparent(0,ii))] +
                       sol [ordUx (quadrant->gparent(1,ii))]);
      Uydof_c  = .5 * (sol [ordUy (quadrant->gparent(0,ii))] +
                       sol [ordUy (quadrant->gparent(1,ii))]);
      Thdof_c  = .5 * (sol [ordTh (quadrant->gparent(0,ii))] +
                       sol [ordTh (quadrant->gparent(1,ii))]);

      Z_node[ii] = .5 * (Z [quadrant->gparent(0,ii)] +
                         Z [quadrant->gparent(1,ii)]);

      P_plus_h_c   = .5 * (P_plus [ordh (quadrant->gparent(0,ii))] +
                           P_plus [ordh (quadrant->gparent(1,ii))]);
      P_minus_h_c  = .5 * (P_minus [ordh (quadrant->gparent(0,ii))] +
                           P_minus [ordh (quadrant->gparent(1,ii))]);

      P_plus_Ux_c   = .5 * (P_plus [ordUx (quadrant->gparent(0,ii))] +
                            P_plus [ordUx (quadrant->gparent(1,ii))]);
      P_minus_Ux_c  = .5 * (P_minus [ordUx (quadrant->gparent(0,ii))] +
                            P_minus [ordUx (quadrant->gparent(1,ii))]);

      P_plus_Uy_c   = .5 * (P_plus [ordUy (quadrant->gparent(0,ii))] +
                            P_plus [ordUy (quadrant->gparent(1,ii))]);
      P_minus_Uy_c  = .5 * (P_minus [ordUy (quadrant->gparent(0,ii))] +
                            P_minus [ordUy (quadrant->gparent(1,ii))]);

      P_plus_Th_c   = .5 * (P_plus [ordTh (quadrant->gparent(0,ii))] +
                            P_plus [ordTh (quadrant->gparent(1,ii))]);
      P_minus_Th_c  = .5 * (P_minus [ordTh (quadrant->gparent(0,ii))] +
                            P_minus [ordTh (quadrant->gparent(1,ii))]);
      
    }

    etadof     [ii] = hdof_c>epsilon ? hdof_c+Z_node[ii]*g_coeff : hdof_c;
    hdof       [ii] = hdof_c;
    Uxdof      [ii] = Uxdof_c;
    Uydof      [ii] = Uydof_c;
    Thdof      [ii] = Thdof_c;

    P_plus_h_dof [ii] = P_plus_h_c;
    P_minus_h_dof[ii] = P_minus_h_c;

    P_plus_Ux_dof [ii] = P_plus_Ux_c;
    P_minus_Ux_dof[ii] = P_minus_Ux_c;

    P_plus_Uy_dof [ii] = P_plus_Uy_c;
    P_minus_Uy_dof[ii] = P_minus_Uy_c;

    P_plus_Th_dof [ii] = P_plus_Th_c;
    P_minus_Th_dof[ii] = P_minus_Th_c;

  }
  
  
  
  // compute local extrema
  const auto h_min_cell  = *std::min_element(etadof.begin(), etadof.end());
  const auto h_max_cell  = *std::max_element(etadof.begin(), etadof.end());
 
  const auto Ux_min_cell = *std::min_element(Uxdof.begin(),  Uxdof.end() );
  const auto Ux_max_cell = *std::max_element(Uxdof.begin(),  Uxdof.end() );

  const auto Uy_min_cell = *std::min_element(Uydof.begin(),  Uydof.end() );
  const auto Uy_max_cell = *std::max_element(Uydof.begin(),  Uydof.end() ); 

  const auto Th_min_cell = *std::min_element(Thdof.begin(),  Thdof.end() );
  const auto Th_max_cell = *std::max_element(Thdof.begin(),  Thdof.end() );  

  bool is_node_in_element = false;

  std::array<double,4> h_min  = {h_min_cell, h_min_cell, h_min_cell, h_min_cell }, h_max  = {h_max_cell, h_max_cell, h_max_cell, h_max_cell },
                       Ux_min = {Ux_min_cell,Ux_min_cell,Ux_min_cell,Ux_min_cell}, Ux_max = {Ux_max_cell,Ux_max_cell,Ux_max_cell,Ux_max_cell},
                       Uy_min = {Uy_min_cell,Uy_min_cell,Uy_min_cell,Uy_min_cell}, Uy_max = {Uy_max_cell,Uy_max_cell,Uy_max_cell,Uy_max_cell},
                       Th_min = {Th_min_cell,Th_min_cell,Th_min_cell,Th_min_cell}, Th_max = {Th_max_cell,Th_max_cell,Th_max_cell,Th_max_cell};


  for (int ii = 0; ii < 4; ++ii){
    for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
         quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
    {

      for (int jj = 0; jj < 4; ++jj) {
        if (xn[ii] == quadrant_nei->p(0, jj) && yn[ii] == quadrant_nei->p(1, jj))
        {
          is_node_in_element = true;
          break;
        }
      }

      if (is_node_in_element)
      {
        for (int jj = 0; jj < 4; ++jj) {

          double eta_current_cell, Z_current_cell, h_current_cell, Ux_current_cell, Uy_current_cell, Th_current_cell;

          if (! quadrant_nei->is_hanging (jj)){
            Z_current_cell  = Z [quadrant_nei->gt (jj)];
            h_current_cell  = sol [ordh  (quadrant_nei->gt (jj))];
            Ux_current_cell = sol [ordUx (quadrant_nei->gt (jj))];
            Uy_current_cell = sol [ordUy (quadrant_nei->gt (jj))];
            Th_current_cell = sol [ordTh (quadrant_nei->gt (jj))];
          } else {
            Z_current_cell  = .5 * (Z [quadrant_nei->gparent(0,jj)] +
                                    Z [quadrant_nei->gparent(1,jj)]);
            h_current_cell  = .5 * (sol [ordh  (quadrant_nei->gparent(0,jj))] +
                                    sol [ordh  (quadrant_nei->gparent(1,jj))]);
            Ux_current_cell = .5 * (sol [ordUx (quadrant_nei->gparent(0,jj))] +
                                    sol [ordUx (quadrant_nei->gparent(1,jj))]);
            Uy_current_cell = .5 * (sol [ordUy (quadrant_nei->gparent(0,jj))] +
                                    sol [ordUy (quadrant_nei->gparent(1,jj))]);
            Th_current_cell = .5 * (sol [ordTh (quadrant_nei->gparent(0,jj))] +
                                    sol [ordTh (quadrant_nei->gparent(1,jj))]);
          }

          eta_current_cell = h_current_cell>epsilon ? h_current_cell+Z_current_cell*g_coeff : h_current_cell;

          h_min[ii]  = std::min(h_min[ii],  eta_current_cell);
          h_max[ii]  = std::max(h_max[ii],  eta_current_cell);

          Ux_min[ii] = std::min(Ux_min[ii], Ux_current_cell);
          Ux_max[ii] = std::max(Ux_max[ii], Ux_current_cell);

          Uy_min[ii] = std::min(Uy_min[ii], Uy_current_cell);
          Uy_max[ii] = std::max(Uy_max[ii], Uy_current_cell);

          Th_min[ii] = std::min(Th_min[ii], Th_current_cell);
          Th_max[ii] = std::max(Th_max[ii], Th_current_cell);

        }

        is_node_in_element = false;
      }


    }
  }



  
  // compute flux correction
  double phi_cell_h = 1., phi_cell_Ux = 1., phi_cell_Uy = 1., phi_cell_Th = 1.;
  for (int ii = 0; ii < 4; ++ii){

    const auto & flux_on_the_node_h  = incr_anti_diff[ordh (index_quadrant)][ii];
    const auto & flux_on_the_node_Ux = incr_anti_diff[ordUx(index_quadrant)][ii];
    const auto & flux_on_the_node_Uy = incr_anti_diff[ordUy(index_quadrant)][ii];
    const auto & flux_on_the_node_Th = incr_anti_diff[ordTh(index_quadrant)][ii];

    const auto hpoint = hdof[ii];

    const auto celerity = std::sqrt(grav*hpoint);
    const auto vel_rusanov_cell_x = hpoint>epsilon ? (std::abs(Uxdof[ii]/hpoint)+celerity) : 0.;
    const auto vel_rusanov_cell_y = hpoint>epsilon ? (std::abs(Uydof[ii]/hpoint)+celerity) : 0.;

    const auto vel_square_rusanov_cell = vel_rusanov_cell_x * vel_rusanov_cell_y;

    flux_limiter(h_min [ii], h_max [ii], etadof[ii], P_plus_h_dof [ii], P_minus_h_dof  [ii], flux_on_the_node_h,  vel_square_rusanov_cell, phi_cell_h );
    flux_limiter(Ux_min[ii], Ux_max[ii], Uxdof [ii], P_plus_Ux_dof[ii], P_minus_Ux_dof [ii], flux_on_the_node_Ux, vel_square_rusanov_cell, phi_cell_Ux);
    flux_limiter(Uy_min[ii], Uy_max[ii], Uydof [ii], P_plus_Uy_dof[ii], P_minus_Uy_dof [ii], flux_on_the_node_Uy, vel_square_rusanov_cell, phi_cell_Uy);
    flux_limiter(Th_min[ii], Th_max[ii], Uydof [ii], P_plus_Th_dof[ii], P_minus_Th_dof [ii], flux_on_the_node_Th, vel_square_rusanov_cell, phi_cell_Th);
  }

  //phi_cell_h = 0., phi_cell_Ux = 0., phi_cell_Uy = 0.;


  for (int ii = 0; ii < 4; ++ii){

    const auto flux_on_the_node_h  = incr_anti_diff[ordh (index_quadrant)][ii]*phi_cell_h;
    const auto flux_on_the_node_Ux = incr_anti_diff[ordUx(index_quadrant)][ii]*phi_cell_Ux;
    const auto flux_on_the_node_Uy = incr_anti_diff[ordUy(index_quadrant)][ii]*phi_cell_Uy;
    const auto flux_on_the_node_Th = incr_anti_diff[ordTh(index_quadrant)][ii]*phi_cell_Th;

    if (! quadrant->is_hanging (ii)){

      incr [ordh  (quadrant->gt (ii))] += flux_on_the_node_h;
      incr [ordUx (quadrant->gt (ii))] += flux_on_the_node_Ux;
      incr [ordUy (quadrant->gt (ii))] += flux_on_the_node_Uy;
      incr [ordTh (quadrant->gt (ii))] += flux_on_the_node_Th;

    } else {

      incr [ordh  (quadrant->gparent(0,ii))] += flux_on_the_node_h; 
      incr [ordh  (quadrant->gparent(1,ii))] += flux_on_the_node_h;
      
      incr [ordUx (quadrant->gparent(0,ii))] += flux_on_the_node_Ux;
      incr [ordUx (quadrant->gparent(1,ii))] += flux_on_the_node_Ux;
      
      incr [ordUy (quadrant->gparent(0,ii))] += flux_on_the_node_Uy;
      incr [ordUy (quadrant->gparent(1,ii))] += flux_on_the_node_Uy;

      incr [ordTh (quadrant->gparent(0,ii))] += flux_on_the_node_Th;
      incr [ordTh (quadrant->gparent(1,ii))] += flux_on_the_node_Th;

    }

  }


}

void
TG2_scheme::flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& vel_square_rusanov_cell, double& phi_cell_Q)
{
  const auto Q_plus  = (Q_max-Q_dof)*dt*vel_square_rusanov_cell;
  const auto Q_minus = (Q_min-Q_dof)*dt*vel_square_rusanov_cell;

  const auto R_plus  = P_plus_Q ==0 ? 1. : std::min(1., Q_plus /P_plus_Q );
  const auto R_minus = P_minus_Q==0 ? 1. : std::min(1., Q_minus/P_minus_Q);

  phi_cell_Q = std::min(phi_cell_Q, flux_on_the_node>=0 ? R_plus : R_minus);
}


void
TG2_scheme::solve_non_lin(const int& kk)
{
  auto & h_c  = sol.get_owned_data ()[kk  ];
  auto & Ux_c = sol.get_owned_data ()[kk+1];
  auto & Uy_c = sol.get_owned_data ()[kk+2];
  auto & Th_c = sol.get_owned_data ()[kk+3];

  const auto & h_c_old  = sold.get_owned_data ()[kk  ];
  const auto & Ux_c_old = sold.get_owned_data ()[kk+1];
  const auto & Uy_c_old = sold.get_owned_data ()[kk+2];
  const auto & Th_c_old = sold.get_owned_data ()[kk+3];


  // solve non-linearities like the first step of the TG2 method to get the complete low order solution,
  h_c += dt*incr.get_owned_data ()[kk]/mass.get_owned_data ()[kk];
  h_c *= (h_c>0.);

  Ux_c += dt*incr.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1];
  Uy_c += dt*incr.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2];

  //Ux_c = (Ux_c + dt*incr.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1] + dt*.5*Ux_src_formula(h_c_old, Ux_c_old, 0., Th_c_old) )/(1.-dt*.5*Ux_src_formula(h_c, 1., 0., Th_c_old));
  //Uy_c = (Uy_c + dt*incr.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] + dt*.5*Uy_src_formula(h_c_old, 0., Uy_c_old, Th_c_old) )/(1.-dt*.5*Uy_src_formula(h_c, 0., 1., Th_c_old));

  Th_c += dt*incr.get_owned_data ()[kk+3]/mass.get_owned_data ()[kk+3];// + dt*.5*Th_src_formula (h_c_old, Ux_c_old, Uy_c_old, Th_c_old);

  //Uy_c = 0.;

  //Newton_energy_balance(h_c, Ux_c, Uy_c, Th_c);

}





void
TG2_scheme::set_dt (const double dt_)
{ dt = dt_; }

void
TG2_scheme::set_old_dt (const double dt_)
{ dt_old = dt_; }

void
TG2_scheme::set_times(const double& t, const double& td, const double& tdd)
{ time = t; timed = td; timedd = tdd; }

double
TG2_scheme::get_dt ()
{ return dt; }


// flux functions
double
TG2_scheme::h_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ 
  // A flux-limiting wetting–drying method for finite-element shallow-water models, with application to the Scheldt Estuary
  return (h>epsilon ? Ux : 0.); 
}

double
TG2_scheme::h_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ 
  // A flux-limiting wetting–drying method for finite-element shallow-water models, with application to the Scheldt Estuary
  return (h>epsilon ? Uy : 0.); 
}

double
TG2_scheme::Ux_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ 
  const auto vel_x = h>epsilon ? Ux/h : 0.;
  return (Ux*vel_x + grav*h*h/2.); 
}
 
double
TG2_scheme::Ux_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return (h>epsilon ? Uy*Ux/h : 0.); }

double
TG2_scheme::Uy_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ return (h>epsilon ? Uy*Ux/h : 0.); }

double
TG2_scheme::Uy_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ 
  const auto vel_y = h>epsilon ? Uy/h : 0.;
  return (Uy*vel_y + grav*h*h/2.); 
}

double
TG2_scheme::Th_flux_formula_x (const double& h, const double& Ux, const double& Uy, const double& Th)
{ return (h>epsilon && !is_isothermal ? Th*Ux/h : 0.); }

double
TG2_scheme::Th_flux_formula_y (const double& h, const double& Ux, const double& Uy, const double& Th)
{ return (h>epsilon && !is_isothermal ? Th*Uy/h : 0.); }



// source terms
double
TG2_scheme::h_src_formula (const double& h, const double& Ux, const double& Uy)
{ 
  return (0.); 
}

double
TG2_scheme::Ux_src_formula (const double& h, const double& Ux, const double& Uy, const double& Th)
{
  const double T = h>epsilon ? Th/h : 0.;
  const double ux = h>epsilon ? Ux/h : 0.;
  const double gamma_fric_over_h = h>epsilon ? 3.*nu_ref/h*std::exp(-b_exp_coeff*(T-T_ref)) : 0.; 

  return ( - gamma_fric_over_h*ux);
}

double
TG2_scheme::Uy_src_formula (const double& h, const double& Ux, const double& Uy, const double& Th)
{
  const double T = h>epsilon ? Th/h : 0.;
  const double uy = h>epsilon ? Uy/h : 0.;
  const double gamma_fric_over_h = h>epsilon ? 3.*nu_ref/h*std::exp(-b_exp_coeff*(T-T_ref)) : 0.; 

  return ( - gamma_fric_over_h*uy);
}



double 
TG2_scheme::src_slope_formula (const double& h, const double& S)
{
  // cell-wise source term to build the incr vector
  return (-grav*S*h);
}


double
TG2_scheme::signum (const double& x)
{ return ((x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0); }


