#include "Taylor_Galerkin.h"
#include <algorithm>


TG2_scheme::TG2_scheme(const Q1& sol,
                       const Q1& sold,
                       const Q1& soldd,
                       Q1& incr,
                       std::vector<std::array<double,4>>& incr_anti_diff,
                       Q1& P_plus,
                       Q1& P_minus,
                       Q0& sol_onehalf,
                       const Q1& mass,
                       const ordering& oh,
                       const ordering& oUx,
                       const ordering& oUy,
                       const Q1& Z,
                       const double& DELTAT,
                       const double& h_min,
                       const bool& is_non_reflBC,
                       const double& density,
                       const double& turbulence_coeff,
                       const double& surface_pressure,
                       const double& bed_friction_angle_rad,
                       const double& fluid_viscosity,
                       const double& yield_shear_stress)
: sol(sol), sold(sold), soldd(soldd), incr(incr), incr_anti_diff(incr_anti_diff), P_plus(P_plus), P_minus(P_minus), sol_onehalf(sol_onehalf), mass(mass),
  ordh(oh), ordUx(oUx), ordUy(oUy), Z(Z), DELTAT(DELTAT), epsilon(h_min), is_non_reflBC(is_non_reflBC), 
  density(density), turbulence_coeff(turbulence_coeff), surface_pressure(surface_pressure), bed_friction_angle_rad(bed_friction_angle_rad), fluid_viscosity(fluid_viscosity), yield_shear_stress(yield_shear_stress)
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
  }
  
  for (int ii = 0; ii < 4; ++ii){
    const auto & hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    
    vel_rusanov_x[ii] = hpoint>epsilon ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
    vel_rusanov_y[ii] = hpoint>epsilon ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;
    
    const auto dtoptx = hpoint>epsilon ? Dx/vel_rusanov_x[ii] : DELTAT;
    const auto dtopty = hpoint>epsilon ? Dy/vel_rusanov_y[ii] : DELTAT;
    const auto dtopt  = dtoptx > dtopty ? dtopty : dtoptx;
    if (dt > dtopt) set_dt (dtopt);
    
//    if (dt<DELTAT)
//    {
//      std::cout << Uxdof[ii] << " " << Uydof[ii] << " " << hpoint << " " << vel_rusanov_x[ii] << " " << vel_rusanov_y[ii] << std::endl;
//      exit( 1.);
//    }
    
//    if (hpoint<0.)
//    {
//      std::cout << hpoint << " " << dt << " " << vel_rusanov_x[ii] << " " << vel_rusanov_y[ii] << std::endl;
//    }
    
    
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
    dh_t = (hdof -  hdof_old)/(time - timed);
    
    h1   = hdof_oldd/((timedd - timed )*(timedd - time ));
    h2   = hdof_old /((timed  - timedd)*(timed  - time ));
    h3   = hdof     /((time   - timedd)*(time   - timed));
    
    a_coeff = h1+h2+h3;
    b_coeff = - (h1*(time+timed) + h2*(time+timedd) + h3*(timed+timedd));
    
    Nu_hmean_cell += std::pow(time-timed,2.)*(4./3.*a_coeff*a_coeff*(time*time+time*timed+timed*timed) + 2.*a_coeff*(b_coeff-dh_t)*(time+timed)+std::pow(b_coeff-dh_t,2.));
    
  }
  Nu_hmean_cell /= 4.;
  nu_htot += Nu_hmean_cell;
  
}

void
TG2_scheme::first_step (tmesh::quadrant_iterator quadrant)
{
  
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii);
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;
  

  double partial_x_Z_average   = 0., partial_y_Z_average    = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      Z_node[ii] = Z[quadrant->gt (ii)];
    }
    else
    {
      Z_node[ii] = .5 * (Z [quadrant->gparent (0, ii)] +
                         Z [quadrant->gparent (1, ii)]);
    }
  }
  partial_x_Z_average = .5*( (Z_node[1]-Z_node[0])/Dx + (Z_node[3]-Z_node[2])/Dx );
  partial_y_Z_average = .5*( (Z_node[2]-Z_node[0])/Dy + (Z_node[3]-Z_node[1])/Dy );


  double h_cell_average        = 0., Ux_cell_average        = 0., Uy_cell_average        = 0.;
  double source_h_cell_average = 0., source_Ux_cell_average = 0., source_Uy_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    double hdof_c, Uxdof_c, Uydof_c;

    if (! quadrant->is_hanging (ii) )
    {
      hdof_c    = sol [ordh  (quadrant->gt (ii) )];
      Uxdof_c   = sol [ordUx (quadrant->gt (ii) )];
      Uydof_c   = sol [ordUy (quadrant->gt (ii) )];
    }
    else
    {
      hdof_c    = .5 * (sol [ordh  (quadrant->gparent (0, ii) )] +
                        sol [ordh  (quadrant->gparent (1, ii) )]);
      Uxdof_c   = .5 * (sol [ordUx (quadrant->gparent (0, ii) )] +
                        sol [ordUx (quadrant->gparent (1, ii) )]);
      Uydof_c   = .5 * (sol [ordUy (quadrant->gparent (0, ii) )] +
                        sol [ordUy (quadrant->gparent (1, ii) )]);
    }

    source_Ux_cell_average += Ux_src_formula (hdof_c, Uxdof_c, Uydof_c, partial_x_Z_average);
    source_Uy_cell_average += Uy_src_formula (hdof_c, Uxdof_c, Uydof_c, partial_y_Z_average);
    
    fluxx_h_node [ii] = h_flux_formula_x   (hdof_c, Uxdof_c, Uydof_c);
    fluxy_h_node [ii] = h_flux_formula_y   (hdof_c, Uxdof_c, Uydof_c);
    fluxx_Ux_node[ii] = Ux_flux_formula_x  (hdof_c, Uxdof_c, Uydof_c);
    fluxy_Ux_node[ii] = Ux_flux_formula_y  (hdof_c, Uxdof_c, Uydof_c);
    fluxx_Uy_node[ii] = Uy_flux_formula_x  (hdof_c, Uxdof_c, Uydof_c);
    fluxy_Uy_node[ii] = Uy_flux_formula_y  (hdof_c, Uxdof_c, Uydof_c);

    h_cell_average  += hdof_c;
    Ux_cell_average += Uxdof_c;
    Uy_cell_average += Uydof_c;

  }
  h_cell_average  /= 4.;
  Ux_cell_average /= 4.;
  Uy_cell_average /= 4.;
  
  source_Ux_cell_average /= 4.;
  source_Uy_cell_average /= 4.;
  
  const auto div_Fh_x = .5*((fluxx_h_node[1]-fluxx_h_node[0]) + (fluxx_h_node[3]-fluxx_h_node[2]));
  const auto div_Fh_y = .5*((fluxy_h_node[2]-fluxy_h_node[0]) + (fluxy_h_node[3]-fluxy_h_node[1]));
  const auto div_Fh_cell = Dy*div_Fh_x + Dx*div_Fh_y;
  
  const auto div_FUx_x = .5*((fluxx_Ux_node[1]-fluxx_Ux_node[0]) + (fluxx_Ux_node[3]-fluxx_Ux_node[2]));
  const auto div_FUx_y = .5*((fluxy_Ux_node[2]-fluxy_Ux_node[0]) + (fluxy_Ux_node[3]-fluxy_Ux_node[1]));
  const auto div_FUx_cell = Dy*div_FUx_x + Dx*div_FUx_y;
  
  const auto div_FUy_x = .5*((fluxx_Uy_node[1]-fluxx_Uy_node[0]) + (fluxx_Uy_node[3]-fluxx_Uy_node[2]));
  const auto div_FUy_y = .5*((fluxy_Uy_node[2]-fluxy_Uy_node[0]) + (fluxy_Uy_node[3]-fluxy_Uy_node[1]));
  const auto div_FUy_cell = Dy*div_FUy_x + Dx*div_FUy_y;
  
  sol_onehalf[ordh  (quadrant->get_forest_quad_idx ())] = h_cell_average  - dt/2.*div_Fh_cell /area;
  sol_onehalf[ordUx (quadrant->get_forest_quad_idx ())] = Ux_cell_average - dt/2.*div_FUx_cell/area + dt/2.*source_Ux_cell_average;
  sol_onehalf[ordUy (quadrant->get_forest_quad_idx ())] = Uy_cell_average - dt/2.*div_FUy_cell/area + dt/2.*source_Uy_cell_average;
  
}

void
TG2_scheme::compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant)
{

  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;


  double h_cell = 0., Ux_cell = 0., Uy_cell = 0.;
  for (int ii = 0; ii < 4; ++ii){

    double hdof_c, Uxdof_c, Uydof_c;

    if (! quadrant->is_hanging (ii)){
      hdof_c    = sol [ordh    (quadrant->gt (ii))];
      Uxdof_c   = sol [ordUx   (quadrant->gt (ii))];
      Uydof_c   = sol [ordUy   (quadrant->gt (ii))];
      
      isdof_or_hanging[ii] = 1.;
    } else {
      hdof_c   = .5 * (sol [ordh  (quadrant->gparent(0,ii))] +
                       sol [ordh  (quadrant->gparent(1,ii))]);
      Uxdof_c  = .5 * (sol [ordUx (quadrant->gparent(0,ii))] +
                       sol [ordUx (quadrant->gparent(1,ii))]);
      Uydof_c  = .5 * (sol [ordUy (quadrant->gparent(0,ii))] +
                       sol [ordUy (quadrant->gparent(1,ii))]);
      
      isdof_or_hanging[ii] = .5;
    }

    hdof [ii] = hdof_c;
    Uxdof[ii] = Uxdof_c;
    Uydof[ii] = Uxdof_c;

    h_cell  += hdof_c;
    Ux_cell += Uxdof_c;
    Uy_cell += Uydof_c;

  }
  h_cell  /= 4.;
  Ux_cell /= 4.;
  Uy_cell /= 4.;


  // weights coefficients for the flux term
  der_coeffs_x = {-Dy/2.*isdof_or_hanging[0], +Dy/2.*isdof_or_hanging[1],
                  -Dy/2.*isdof_or_hanging[2], +Dy/2.*isdof_or_hanging[3]};
  
  der_coeffs_y = {-Dx/2.*isdof_or_hanging[0], -Dx/2.*isdof_or_hanging[1],
                  +Dx/2.*isdof_or_hanging[2], +Dx/2.*isdof_or_hanging[3]};

  double vel_rusanov_cell_x = 0., vel_rusanov_cell_y = 0.;
  for (int ii = 0; ii < 4; ++ii){
    const auto & hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    vel_rusanov_cell_x += hpoint>epsilon ? (std::abs(Uxdof[ii]/hpoint)+celerity) : 0.;
    vel_rusanov_cell_y += hpoint>epsilon ? (std::abs(Uydof[ii]/hpoint)+celerity) : 0.;
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.;
  

  
  grad_cell_h  = {.5 * ( (hdof [3] - hdof [2]) + (hdof [1] - hdof [0]) ), .5 * ( (hdof [2] - hdof [0]) + (hdof [3] - hdof [1]) )};
  grad_cell_Ux = {.5 * ( (Uxdof[3] - Uxdof[2]) + (Uxdof[1] - Uxdof[0]) ), .5 * ( (Uxdof[2] - Uxdof[0]) + (Uxdof[3] - Uxdof[1]) )};
  grad_cell_Uy = {.5 * ( (Uydof[3] - Uydof[2]) + (Uydof[1] - Uydof[0]) ), .5 * ( (Uydof[2] - Uydof[0]) + (Uydof[3] - Uydof[1]) )};


  const auto & h_cell_12  = sol_onehalf[ordh  (quadrant->get_forest_quad_idx ())]; 
  const auto & Ux_cell_12 = sol_onehalf[ordUx (quadrant->get_forest_quad_idx ())]; 
  const auto & Uy_cell_12 = sol_onehalf[ordUy (quadrant->get_forest_quad_idx ())]; 


  const auto diff_term_h_x  = .5 * grad_cell_h [0] * vel_rusanov_cell_y;
  const auto diff_term_h_y  = .5 * grad_cell_h [1] * vel_rusanov_cell_x;

  const auto diff_term_Ux_x = .5 * grad_cell_Ux[0] * vel_rusanov_cell_y;
  const auto diff_term_Ux_y = .5 * grad_cell_Ux[1] * vel_rusanov_cell_x;

  const auto diff_term_Uy_x = .5 * grad_cell_Uy[0] * vel_rusanov_cell_y;
  const auto diff_term_Uy_y = .5 * grad_cell_Uy[1] * vel_rusanov_cell_x;

  // low order flux
  const auto F_star_h_x_l  = h_flux_formula_x (h_cell, Ux_cell, Uy_cell) - diff_term_h_x;
  const auto F_star_h_y_l  = h_flux_formula_y (h_cell, Ux_cell, Uy_cell) - diff_term_h_y;

  const auto F_star_Ux_x_l = Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_x;
  const auto F_star_Ux_y_l = Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_y;

  const auto F_star_Uy_x_l = Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_x;
  const auto F_star_Uy_y_l = Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_y;

  // antidiffusive nodal flux
  const auto F_star_h_x  = h_flux_formula_x (h_cell_12, Ux_cell_12, Uy_cell_12) - F_star_h_x_l;
  const auto F_star_h_y  = h_flux_formula_y (h_cell_12, Ux_cell_12, Uy_cell_12) - F_star_h_y_l;

  const auto F_star_Ux_x = Ux_flux_formula_x(h_cell_12, Ux_cell_12, Uy_cell_12) - F_star_Ux_x_l;
  const auto F_star_Ux_y = Ux_flux_formula_y(h_cell_12, Ux_cell_12, Uy_cell_12) - F_star_Ux_y_l;

  const auto F_star_Uy_x = Uy_flux_formula_x(h_cell_12, Ux_cell_12, Uy_cell_12) - F_star_Uy_x_l;
  const auto F_star_Uy_y = Uy_flux_formula_y(h_cell_12, Ux_cell_12, Uy_cell_12) - F_star_Uy_y_l;



  double partial_x_Z_average = 0., partial_y_Z_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {

    if (! quadrant->is_hanging (ii) )
    {
      Z_node[ii] = Z[quadrant->gt (ii)];
    }
    else
    {
      Z_node[ii] = .5 * (Z [quadrant->gparent (0, ii)] +
                         Z [quadrant->gparent (1, ii)]);
    }
  }
  partial_x_Z_average = .5*( (Z_node[1]-Z_node[0])/Dx + (Z_node[3]-Z_node[2])/Dx );
  partial_y_Z_average = .5*( (Z_node[2]-Z_node[0])/Dy + (Z_node[3]-Z_node[1])/Dy );



  for (int ii = 0; ii < 4; ++ii){

    const auto h_  = der_coeffs_x[ii] * F_star_h_x  + der_coeffs_y[ii] * F_star_h_y;
    const auto Ux_ = der_coeffs_x[ii] * F_star_Ux_x + der_coeffs_y[ii] * F_star_Ux_y;
    const auto Uy_ = der_coeffs_x[ii] * F_star_Uy_x + der_coeffs_y[ii] * F_star_Uy_y;

    const auto h_incr  = der_coeffs_x[ii] * F_star_h_x_l  + der_coeffs_y[ii] * F_star_h_y_l;
    const auto Ux_incr = der_coeffs_x[ii] * F_star_Ux_x_l + der_coeffs_y[ii] * F_star_Ux_y_l + .25*area*isdof_or_hanging[ii] * Ux_src_formula(h_cell_12, Ux_cell_12, Uy_cell_12, partial_x_Z_average);
    const auto Uy_incr = der_coeffs_x[ii] * F_star_Uy_x_l + der_coeffs_y[ii] * F_star_Uy_y_l + .25*area*isdof_or_hanging[ii] * Uy_src_formula(h_cell_12, Ux_cell_12, Uy_cell_12, partial_y_Z_average);
    
    incr_anti_diff[ordh (quadrant->get_forest_quad_idx ())][ii] = h_;
    incr_anti_diff[ordUx(quadrant->get_forest_quad_idx ())][ii] = Ux_;
    incr_anti_diff[ordUy(quadrant->get_forest_quad_idx ())][ii] = Uy_;

    //std::cout << h_incr << " " << Ux_incr << " " << Uy_incr << " " << h_cell << " " << Ux_cell << " " << Uy_cell << " " << diff_term_h_x << std::endl;

    if (! quadrant->is_hanging (ii)){

      incr [ordh  (quadrant->gt (ii))] += h_incr;
      incr [ordUx (quadrant->gt (ii))] += Ux_incr;
      incr [ordUy (quadrant->gt (ii))] += Uy_incr;


      P_plus [ordh  (quadrant->gt (ii))] += std::max(0., h_ );
      P_plus [ordUx (quadrant->gt (ii))] += std::max(0., Ux_);
      P_plus [ordUy (quadrant->gt (ii))] += std::max(0., Uy_);

      P_minus [ordh  (quadrant->gt (ii))] += std::min(0., h_ );
      P_minus [ordUx (quadrant->gt (ii))] += std::min(0., Ux_);
      P_minus [ordUy (quadrant->gt (ii))] += std::min(0., Uy_);
      
      
      const auto boundary_idxx = quadrant->ex (ii);
      const auto boundary_idxy = quadrant->ey (ii);
      
      // if (boundary_idxx == tmesh::quadrant_t::NOT_ON_BOUNDARY && boundary_idxy == tmesh::quadrant_t::NOT_ON_BOUNDARY && quadrant->gt (ii)==424)
      // {
      //   std::cout << h_incr << " " << Ux_incr << " " << Uy_incr << " " << quadrant->gt (ii) << std::endl;
      // }

      // if (boundary_idxy == tmesh::quadrant_t::NOT_ON_BOUNDARY && quadrant->gt (ii)==560)
      // {
      //   std::cout << h_incr << " " << Ux_incr << " " << Uy_incr << " " << quadrant->gt (ii) << std::endl;
      // }


      double F_star_h_x_b_l  = 0., F_star_h_y_b_l  = 0.,
             F_star_Ux_x_b_l = 0., F_star_Ux_y_b_l = 0., 
             F_star_Uy_x_b_l = 0., F_star_Uy_y_b_l = 0.;

      double F_star_h_x_b  = 0., F_star_h_y_b  = 0.,
             F_star_Ux_x_b = 0., F_star_Ux_y_b = 0., 
             F_star_Uy_x_b = 0., F_star_Uy_y_b = 0.;


      if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY) {

        F_star_h_x_b_l  = h_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_h_x;
        F_star_h_y_b_l  = h_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_h_y;

        F_star_Ux_x_b_l = Ux_flux_formula_x(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_Ux_x;
        F_star_Ux_y_b_l = Ux_flux_formula_y(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_Ux_y;

        F_star_Uy_x_b_l = Uy_flux_formula_x(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_Uy_x;
        F_star_Uy_y_b_l = Uy_flux_formula_y(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_Uy_y;
    

        F_star_h_x_b  = h_flux_formula_x (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_h_x_b_l;
        F_star_h_y_b  = h_flux_formula_y (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_h_y_b_l;

        F_star_Ux_x_b = Ux_flux_formula_x(h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_Ux_x_b_l;
        F_star_Ux_y_b = Ux_flux_formula_y(h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_Ux_y_b_l;

        F_star_Uy_x_b = Uy_flux_formula_x(h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_Uy_x_b_l;
        F_star_Uy_y_b = Uy_flux_formula_y(h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_Uy_y_b_l;

        
        auto h_b_l  = -der_coeffs_x[ii] * F_star_h_x_b_l  + der_coeffs_y[ii] * F_star_h_y_b_l;
        auto Ux_b_l = -der_coeffs_x[ii] * F_star_Ux_x_b_l + der_coeffs_y[ii] * F_star_Ux_y_b_l;
        auto Uy_b_l = -der_coeffs_x[ii] * F_star_Uy_x_b_l + der_coeffs_y[ii] * F_star_Uy_y_b_l;

        auto h_b  = -der_coeffs_x[ii] * F_star_h_x_b  + der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = -der_coeffs_x[ii] * F_star_Ux_x_b + der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = -der_coeffs_x[ii] * F_star_Uy_x_b + der_coeffs_y[ii] * F_star_Uy_y_b;
        
        

        incr [ordh  (quadrant->gt (ii))] += h_b_l;
        incr [ordUx (quadrant->gt (ii))] += Ux_b_l;
        incr [ordUy (quadrant->gt (ii))] += Uy_b_l;


        incr_anti_diff[ordh (quadrant->get_forest_quad_idx ())][ii] += h_b;
        incr_anti_diff[ordUx(quadrant->get_forest_quad_idx ())][ii] += Ux_b;
        incr_anti_diff[ordUy(quadrant->get_forest_quad_idx ())][ii] += Uy_b;

        P_plus [ordh    (quadrant->gt (ii))] += std::max(0., h_b );
        P_plus [ordUx   (quadrant->gt (ii))] += std::max(0., Ux_b);
        P_plus [ordUy   (quadrant->gt (ii))] += std::max(0., Uy_b);

        P_minus [ordh    (quadrant->gt (ii))] += std::min(0., h_b );
        P_minus [ordUx   (quadrant->gt (ii))] += std::min(0., Ux_b);
        P_minus [ordUy   (quadrant->gt (ii))] += std::min(0., Uy_b);
        
        
        for (int jj = 0; jj < 4; ++jj){
          if ( quadrant->is_hanging (jj) && (der_coeffs_y[jj]*der_coeffs_y[ii] > 0.) )
          {
            F_star_h_y_b_l  = h_flux_formula_y  (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_h_y;
            F_star_Ux_y_b_l = Ux_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_Ux_y;
            F_star_Uy_y_b_l = Uy_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell) - diff_term_Uy_y;

            F_star_h_y_b  = h_flux_formula_y (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_h_y_b_l;
            F_star_Ux_y_b = Ux_flux_formula_y(h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_Ux_y_b_l;
            F_star_Uy_y_b = Uy_flux_formula_y(h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, Uy_cell_12) - F_star_Uy_y_b_l;

            
            h_b_l  = 2. * der_coeffs_y[jj] * F_star_h_y_b_l;
            Ux_b_l = 2. * der_coeffs_y[jj] * F_star_Ux_y_b_l;  
            Uy_b_l = 2. * der_coeffs_y[jj] * F_star_Uy_y_b_l;

            h_b  = 2. * der_coeffs_y[jj] * F_star_h_y_b;
            Ux_b = 2. * der_coeffs_y[jj] * F_star_Ux_y_b; 
            Uy_b = 2. * der_coeffs_y[jj] * F_star_Uy_y_b;


            incr [ordh  (quadrant->gt (ii))] += h_b_l;
            incr [ordUx (quadrant->gt (ii))] += Ux_b_l;
            incr [ordUy (quadrant->gt (ii))] += Uy_b_l;

            incr_anti_diff[ordh (quadrant->get_forest_quad_idx ())][ii] += h_b;
            incr_anti_diff[ordUx(quadrant->get_forest_quad_idx ())][ii] += Ux_b;
            incr_anti_diff[ordUy(quadrant->get_forest_quad_idx ())][ii] += Uy_b;
            
            P_plus [ordh  (quadrant->gt (ii))] += std::max(0., h_b );
            P_plus [ordUx (quadrant->gt (ii))] += std::max(0., Ux_b); 
            P_plus [ordUy (quadrant->gt (ii))] += std::max(0., Uy_b);

            P_minus [ordh  (quadrant->gt (ii))] += std::min(0., h_b );
            P_minus [ordUx (quadrant->gt (ii))] += std::min(0., Ux_b);
            P_minus [ordUy (quadrant->gt (ii))] += std::min(0., Uy_b);
          }
          
        }
        
        
        
        
      }
      
      if (boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY) {

        F_star_h_x_b_l  = h_flux_formula_x (h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_h_x;
        F_star_h_y_b_l  = h_flux_formula_y (h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_h_y;

        F_star_Ux_x_b_l = Ux_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Ux_x;
        F_star_Ux_y_b_l = Ux_flux_formula_y(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Ux_y;

        F_star_Uy_x_b_l = Uy_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Uy_x;
        F_star_Uy_y_b_l = Uy_flux_formula_y(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Uy_y;


        F_star_h_x_b  = h_flux_formula_x (h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_h_x_b_l;
        F_star_h_y_b  = h_flux_formula_y (h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_h_y_b_l;

        F_star_Ux_x_b = Ux_flux_formula_x(h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Ux_x_b_l;
        F_star_Ux_y_b = Ux_flux_formula_y(h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Ux_y_b_l;

        F_star_Uy_x_b = Uy_flux_formula_x(h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Uy_x_b_l;
        F_star_Uy_y_b = Uy_flux_formula_y(h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Uy_y_b_l;

        
        auto h_b_l  = der_coeffs_x[ii] * F_star_h_x_b_l  - der_coeffs_y[ii] * F_star_h_y_b_l;
        auto Ux_b_l = der_coeffs_x[ii] * F_star_Ux_x_b_l - der_coeffs_y[ii] * F_star_Ux_y_b_l;
        auto Uy_b_l = der_coeffs_x[ii] * F_star_Uy_x_b_l - der_coeffs_y[ii] * F_star_Uy_y_b_l;

        auto h_b  = der_coeffs_x[ii] * F_star_h_x_b  - der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = der_coeffs_x[ii] * F_star_Ux_x_b - der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = der_coeffs_x[ii] * F_star_Uy_x_b - der_coeffs_y[ii] * F_star_Uy_y_b;
        
        //std::cout << h_incr << " " << Ux_incr << " " << Uy_incr << " " << h_b_l << " " << Ux_b_l << " " << Uy_b_l << std::endl;

        incr [ordh  (quadrant->gt (ii))] += h_b_l;
        incr [ordUx (quadrant->gt (ii))] += Ux_b_l;
        incr [ordUy (quadrant->gt (ii))] += Uy_b_l;

        incr_anti_diff[ordh (quadrant->get_forest_quad_idx ())][ii] += h_b;
        incr_anti_diff[ordUx(quadrant->get_forest_quad_idx ())][ii] += Ux_b;
        incr_anti_diff[ordUy(quadrant->get_forest_quad_idx ())][ii] += Uy_b;

        P_plus [ordh    (quadrant->gt (ii))] += std::max(0., h_b );
        P_plus [ordUx   (quadrant->gt (ii))] += std::max(0., Ux_b);
        P_plus [ordUy   (quadrant->gt (ii))] += std::max(0., Uy_b);

        P_minus [ordh    (quadrant->gt (ii))] += std::min(0., h_b );
        P_minus [ordUx   (quadrant->gt (ii))] += std::min(0., Ux_b);
        P_minus [ordUy   (quadrant->gt (ii))] += std::min(0., Uy_b);
        
        
        for (int jj = 0; jj < 4; ++jj){
          if ( quadrant->is_hanging (jj) && (der_coeffs_x[jj]*der_coeffs_x[ii] > 0.) )
          {
            F_star_h_x_b_l  = h_flux_formula_x (h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_h_x;
            F_star_Ux_x_b_l = Ux_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Ux_x;
            F_star_Uy_x_b_l = Uy_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Uy_x;

            F_star_h_x_b  = h_flux_formula_x (h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_h_x_b_l;
            F_star_Ux_x_b = Ux_flux_formula_x(h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Ux_x_b_l;
            F_star_Uy_x_b = Uy_flux_formula_x(h_cell_12, Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Uy_x_b_l;
            
            h_b_l  = 2. * der_coeffs_x[jj] * F_star_h_x_b_l;
            Ux_b_l = 2. * der_coeffs_x[jj] * F_star_Ux_x_b_l;
            Uy_b_l = 2. * der_coeffs_x[jj] * F_star_Uy_x_b_l;

            h_b  = 2. * der_coeffs_x[jj] * F_star_h_x_b;
            Ux_b = 2. * der_coeffs_x[jj] * F_star_Ux_x_b;
            Uy_b = 2. * der_coeffs_x[jj] * F_star_Uy_x_b;
            

            incr [ordh  (quadrant->gt (ii))] += h_b_l;
            incr [ordUx (quadrant->gt (ii))] += Ux_b_l;
            incr [ordUy (quadrant->gt (ii))] += Uy_b_l;


            incr_anti_diff[ordh (quadrant->get_forest_quad_idx ())][ii] += h_b;
            incr_anti_diff[ordUx(quadrant->get_forest_quad_idx ())][ii] += Ux_b;
            incr_anti_diff[ordUy(quadrant->get_forest_quad_idx ())][ii] += Uy_b;

            P_plus [ordh    (quadrant->gt (ii))] += std::max(0., h_b );
            P_plus [ordUx   (quadrant->gt (ii))] += std::max(0., Ux_b);
            P_plus [ordUy   (quadrant->gt (ii))] += std::max(0., Uy_b);

            P_minus [ordh    (quadrant->gt (ii))] += std::min(0., h_b );
            P_minus [ordUx   (quadrant->gt (ii))] += std::min(0., Ux_b);
            P_minus [ordUy   (quadrant->gt (ii))] += std::min(0., Uy_b);
          }
          
        }
        
        
      }
      
      // corner points!
      if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY && boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY)
      {
        
        F_star_h_x_b_l  = h_flux_formula_x  (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_h_x;
        F_star_h_y_b_l  = h_flux_formula_y  (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_h_y;
        
        F_star_Ux_x_b_l = Ux_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Ux_x;
        F_star_Ux_y_b_l = Ux_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Ux_y;
        
        F_star_Uy_x_b_l = Uy_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Uy_x;
        F_star_Uy_y_b_l = Uy_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell) - diff_term_Uy_y;


        F_star_h_x_b  = h_flux_formula_x  (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_h_x_b_l;
        F_star_h_y_b  = h_flux_formula_y  (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_h_y_b_l;
        
        F_star_Ux_x_b = Ux_flux_formula_x  (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Ux_x_b_l;
        F_star_Ux_y_b = Ux_flux_formula_y  (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Ux_y_b_l;
        
        F_star_Uy_x_b = Uy_flux_formula_x  (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Uy_x_b_l;
        F_star_Uy_y_b = Uy_flux_formula_y  (h_cell_12, is_non_reflBC ? Ux_cell_12 : -Ux_cell_12, is_non_reflBC ? Uy_cell_12 : -Uy_cell_12) - F_star_Uy_y_b_l;
        
        
        auto h_b_l  = -der_coeffs_x[ii] * F_star_h_x_b_l  - der_coeffs_y[ii] * F_star_h_y_b_l;
        auto Ux_b_l = -der_coeffs_x[ii] * F_star_Ux_x_b_l - der_coeffs_y[ii] * F_star_Ux_y_b_l;
        auto Uy_b_l = -der_coeffs_x[ii] * F_star_Uy_x_b_l - der_coeffs_y[ii] * F_star_Uy_y_b_l;

        auto h_b  = -der_coeffs_x[ii] * F_star_h_x_b  - der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = -der_coeffs_x[ii] * F_star_Ux_x_b - der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = -der_coeffs_x[ii] * F_star_Uy_x_b - der_coeffs_y[ii] * F_star_Uy_y_b;
        

        incr [ordh  (quadrant->gt (ii))] += h_b_l;
        incr [ordUx (quadrant->gt (ii))] += Ux_b_l;
        incr [ordUy (quadrant->gt (ii))] += Uy_b_l;


        incr_anti_diff[ordh (quadrant->get_forest_quad_idx ())][ii] += h_b;
        incr_anti_diff[ordUx(quadrant->get_forest_quad_idx ())][ii] += Ux_b;
        incr_anti_diff[ordUy(quadrant->get_forest_quad_idx ())][ii] += Uy_b;

        P_plus [ordh    (quadrant->gt (ii))] += std::max(0., h_b );
        P_plus [ordUx   (quadrant->gt (ii))] += std::max(0., Ux_b);
        P_plus [ordUy   (quadrant->gt (ii))] += std::max(0., Uy_b);

        P_minus [ordh    (quadrant->gt (ii))] += std::min(0., h_b );
        P_minus [ordUx   (quadrant->gt (ii))] += std::min(0., Ux_b);
        P_minus [ordUy   (quadrant->gt (ii))] += std::min(0., Uy_b);
        
      }
      
      
    } else {

      incr [ordh  (quadrant->gparent(0,ii))] += h_incr;
      incr [ordh  (quadrant->gparent(1,ii))] += h_incr;
      
      incr [ordUx (quadrant->gparent(0,ii))] += Ux_incr;
      incr [ordUx (quadrant->gparent(1,ii))] += Ux_incr;
      
      incr [ordUy (quadrant->gparent(0,ii))] += Uy_incr;
      incr [ordUy (quadrant->gparent(1,ii))] += Uy_incr;


      
      P_plus [ordh  (quadrant->gparent(0,ii))] += std::max(0., h_);
      P_plus [ordh  (quadrant->gparent(1,ii))] += std::max(0., h_);
      
      P_plus [ordUx (quadrant->gparent(0,ii))] += std::max(0., Ux_);
      P_plus [ordUx (quadrant->gparent(1,ii))] += std::max(0., Ux_);
      
      P_plus [ordUy (quadrant->gparent(0,ii))] += std::max(0., Uy_);
      P_plus [ordUy (quadrant->gparent(1,ii))] += std::max(0., Uy_);



      P_minus [ordh  (quadrant->gparent(0,ii))] += std::min(0., h_);
      P_minus [ordh  (quadrant->gparent(1,ii))] += std::min(0., h_);
      
      P_minus [ordUx (quadrant->gparent(0,ii))] += std::min(0., Ux_);
      P_minus [ordUx (quadrant->gparent(1,ii))] += std::min(0., Ux_);
      
      P_minus [ordUy (quadrant->gparent(0,ii))] += std::min(0., Uy_);
      P_minus [ordUy (quadrant->gparent(1,ii))] += std::min(0., Uy_);
      
    }
  }

}





void
TG2_scheme::second_step (tmesh::quadrant_iterator quadrant)
{
  
  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  
  for (int ii = 0; ii < 4; ++ii){

    double hdof_c, Uxdof_c, Uydof_c, P_plus_h_c, P_minus_h_c, P_plus_Ux_c, P_minus_Ux_c, P_plus_Uy_c, P_minus_Uy_c;

    if (! quadrant->is_hanging (ii)){
      hdof_c      = sol [ordh    (quadrant->gt (ii))];
      Uxdof_c     = sol [ordUx   (quadrant->gt (ii))];
      Uydof_c     = sol [ordUy   (quadrant->gt (ii))];

      P_plus_h_c   = P_plus [ordh   (quadrant->gt (ii))];
      P_minus_h_c  = P_minus[ordh   (quadrant->gt (ii))];

      P_plus_Ux_c   = P_plus [ordUx   (quadrant->gt (ii))];
      P_minus_Ux_c  = P_minus[ordUx   (quadrant->gt (ii))];

      P_plus_Uy_c   = P_plus [ordUy   (quadrant->gt (ii))];
      P_minus_Uy_c  = P_minus[ordUy   (quadrant->gt (ii))];
      

    } else {
      hdof_c   = .5 * (sol [ordh  (quadrant->gparent(0,ii))] +
                       sol [ordh  (quadrant->gparent(1,ii))]);
      Uxdof_c  = .5 * (sol [ordUx (quadrant->gparent(0,ii))] +
                       sol [ordUx (quadrant->gparent(1,ii))]);
      Uydof_c  = .5 * (sol [ordUy (quadrant->gparent(0,ii))] +
                       sol [ordUy (quadrant->gparent(1,ii))]);

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
      
    }

    hdof       [ii] = hdof_c;
    Uxdof      [ii] = Uxdof_c;
    Uydof      [ii] = Uxdof_c;

    P_plus_h_dof [ii] = P_plus_h_c;
    P_minus_h_dof[ii] = P_minus_h_c;

    P_plus_Ux_dof [ii] = P_plus_Ux_c;
    P_minus_Ux_dof[ii] = P_minus_Ux_c;

    P_plus_Uy_dof [ii] = P_plus_Uy_c;
    P_minus_Uy_dof[ii] = P_minus_Uy_c;

  }
  
  
  
  // compute local extrema
  const auto h_min_cell  = *std::min_element(hdof.begin(),  hdof.end() );
  const auto h_max_cell  = *std::max_element(hdof.begin(),  hdof.end() );
 
  const auto Ux_min_cell = *std::min_element(Uxdof.begin(), Uxdof.end());
  const auto Ux_max_cell = *std::max_element(Uxdof.begin(), Uxdof.end());

  const auto Uy_min_cell = *std::min_element(Uydof.begin(), Uydof.end());
  const auto Uy_max_cell = *std::max_element(Uydof.begin(), Uydof.end());  

  bool is_node_in_element = false;

  std::array<double,4> h_min  = {h_min_cell, h_min_cell, h_min_cell, h_min_cell }, h_max  = {h_max_cell, h_max_cell, h_max_cell, h_max_cell },
                       Ux_min = {Ux_min_cell,Ux_min_cell,Ux_min_cell,Ux_min_cell}, Ux_max = {Ux_max_cell,Ux_max_cell,Ux_max_cell,Ux_max_cell},
                       Uy_min = {Uy_min_cell,Uy_min_cell,Uy_min_cell,Uy_min_cell}, Uy_max = {Uy_max_cell,Uy_max_cell,Uy_max_cell,Uy_max_cell};

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

          double h_current_cell, Ux_current_cell, Uy_current_cell;

          if (! quadrant_nei->is_hanging (jj)){
            h_current_cell  = sol [ordh  (quadrant_nei->gt (jj))];
            Ux_current_cell = sol [ordUx (quadrant_nei->gt (jj))];
            Uy_current_cell = sol [ordUy (quadrant_nei->gt (jj))];
          } else {
            h_current_cell  = .5 * (sol [ordh  (quadrant_nei->gparent(0,jj))] +
                                    sol [ordh  (quadrant_nei->gparent(1,jj))]);
            Ux_current_cell = .5 * (sol [ordUx (quadrant_nei->gparent(0,jj))] +
                                    sol [ordUx (quadrant_nei->gparent(1,jj))]);
            Uy_current_cell = .5 * (sol [ordUy (quadrant_nei->gparent(0,jj))] +
                                    sol [ordUy (quadrant_nei->gparent(1,jj))]);
          }

          h_min[ii]  = std::min(h_min[ii],  h_current_cell);
          h_max[ii]  = std::max(h_max[ii],  h_current_cell);

          Ux_min[ii] = std::min(Ux_min[ii], Ux_current_cell);
          Ux_max[ii] = std::max(Ux_max[ii], Ux_current_cell);

          Uy_min[ii] = std::min(Uy_min[ii], Uy_current_cell);
          Uy_max[ii] = std::max(Uy_max[ii], Uy_current_cell);

        }

        is_node_in_element = false;
      }


    }
  }



  
  // compute flux limiter, grad limiter
  double phi_cell_h = 1., phi_cell_Ux = 1., phi_cell_Uy = 1.;
  for (int ii = 0; ii < 4; ++ii){

    const auto & flux_on_the_node_h  = incr_anti_diff[ordh (quadrant->get_forest_quad_idx ())][ii];
    const auto & flux_on_the_node_Ux = incr_anti_diff[ordUx(quadrant->get_forest_quad_idx ())][ii];
    const auto & flux_on_the_node_Uy = incr_anti_diff[ordUy(quadrant->get_forest_quad_idx ())][ii];


    double mass_node_h, mass_node_Ux, mass_node_Uy;
    if (! quadrant->is_hanging (ii)){
      mass_node_h  = mass[ordh (quadrant->gt(ii))];
      mass_node_Ux = mass[ordUx(quadrant->gt(ii))];
      mass_node_Uy = mass[ordUy(quadrant->gt(ii))];
    } else {
      mass_node_h  = .5 * (mass [ordh   (quadrant->gparent(0,ii))] +
                           mass [ordh   (quadrant->gparent(1,ii))]);
      mass_node_Ux = .5 * (mass [ordUx  (quadrant->gparent(0,ii))] +
                           mass [ordUx  (quadrant->gparent(1,ii))]);
      mass_node_Uy = .5 * (mass [ordUy  (quadrant->gparent(0,ii))] +
                           mass [ordUy  (quadrant->gparent(1,ii))]);
    }


    const auto & hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    const auto vel_rusanov_cell_x = hpoint>epsilon ? (std::abs(Uxdof[ii]/hpoint)+celerity) : 0.;
    const auto vel_rusanov_cell_y = hpoint>epsilon ? (std::abs(Uydof[ii]/hpoint)+celerity) : 0.;

    const auto vel_rusanov_cell = vel_rusanov_cell_x * vel_rusanov_cell_y;

    flux_limiter(h_min [ii], h_max [ii], hdof [ii], P_plus_h_dof [ii], P_minus_h_dof  [ii], flux_on_the_node_h,  vel_rusanov_cell,  phi_cell_h );
    flux_limiter(Ux_min[ii], Ux_max[ii], Uxdof[ii], P_plus_Ux_dof[ii], P_minus_Ux_dof [ii], flux_on_the_node_Ux, vel_rusanov_cell, phi_cell_Ux);
    flux_limiter(Uy_min[ii], Uy_max[ii], Uydof[ii], P_plus_Uy_dof[ii], P_minus_Uy_dof [ii], flux_on_the_node_Uy, vel_rusanov_cell, phi_cell_Uy);
  }

  //std::cout << phi_cell_h << " " << phi_cell_Ux << " " << phi_cell_Uy << std::endl;
  //phi_cell_h = 0.; phi_cell_Ux = 0.; phi_cell_Uy = 0.;

  for (int ii = 0; ii < 4; ++ii){
    const auto flux_on_the_node_h  = incr_anti_diff[ordh (quadrant->get_forest_quad_idx ())][ii]*phi_cell_h;
    const auto flux_on_the_node_Ux = incr_anti_diff[ordUx(quadrant->get_forest_quad_idx ())][ii]*phi_cell_Ux;
    const auto flux_on_the_node_Uy = incr_anti_diff[ordUy(quadrant->get_forest_quad_idx ())][ii]*phi_cell_Uy;

    if (! quadrant->is_hanging (ii)){

      incr [ordh  (quadrant->gt (ii))] += flux_on_the_node_h;
      incr [ordUx (quadrant->gt (ii))] += flux_on_the_node_Ux;
      incr [ordUy (quadrant->gt (ii))] += flux_on_the_node_Uy;

    } else {

      incr [ordh  (quadrant->gparent(0,ii))] += flux_on_the_node_h;
      incr [ordh  (quadrant->gparent(1,ii))] += flux_on_the_node_h;
      
      incr [ordUx (quadrant->gparent(0,ii))] += flux_on_the_node_Ux;
      incr [ordUx (quadrant->gparent(1,ii))] += flux_on_the_node_Ux;
      
      incr [ordUy (quadrant->gparent(0,ii))] += flux_on_the_node_Uy;
      incr [ordUy (quadrant->gparent(1,ii))] += flux_on_the_node_Uy;

    }

  }



}

void
TG2_scheme::flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& vel_rusanov_cell, double& phi_cell_Q)
{
  const auto Q_plus  = (Q_max-Q_dof)*dt*vel_rusanov_cell;
  const auto Q_minus = (Q_min-Q_dof)*dt*vel_rusanov_cell;

  const auto R_plus  = P_plus_Q ==0 ? 1 : std::min(1., Q_plus /P_plus_Q );
  const auto R_minus = P_minus_Q==0 ? 1 : std::min(1., Q_minus/P_minus_Q);

  phi_cell_Q = std::min(phi_cell_Q, flux_on_the_node>=0 ? R_plus : R_minus);

}



void
TG2_scheme::set_dt (const double dt_)
{ dt = dt_; }

void
TG2_scheme::set_times(const double& t, const double& td, const double& tdd)
{ time = t; timed = td; timedd = tdd; }

double
TG2_scheme::get_dt ()
{ return dt; }


// flux functions
double
TG2_scheme::h_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ? Ux : 0.; }

double
TG2_scheme::h_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ? Uy : 0.; }

double
TG2_scheme::Ux_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ?  (Ux*Ux/h + grav*h*h/2.) : 0.; }

double
TG2_scheme::Ux_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ? (Uy*Ux/h) : 0.; }

double
TG2_scheme::Uy_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ? (Uy*Ux/h) : 0.; }

double
TG2_scheme::Uy_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ? (Uy*Uy/h + grav*h*h/2.) : 0.; }


// source terms
double
TG2_scheme::h_src_formula (const double& h, const double& Ux, const double& Uy)
{ return (0.); }

double
TG2_scheme::Ux_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdx)
{
  const double bed_pressure = density*grav*h - surface_pressure;
  const double abs_vel = h > epsilon ? std::sqrt( std::pow((Ux/h),2.) + std::pow((Uy/h),2.) ) : 0.;
  return (-grav*h*dZdx - (grav*abs_vel/turbulence_coeff + h > epsilon ? bed_pressure*std::tan(bed_friction_angle_rad)/abs_vel/density : 0.)*(h>epsilon ? Ux/h : 0.) );
}

double
TG2_scheme::Uy_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdy)
{
  const double bed_pressure = density*grav*h - surface_pressure;
  const double abs_vel = h > epsilon ? std::sqrt( std::pow((Ux/h),2.) + std::pow((Uy/h),2.) ) : 0.;
  return (-grav*h*dZdy - (grav*abs_vel/turbulence_coeff + h > epsilon ? bed_pressure*std::tan(bed_friction_angle_rad)/abs_vel/density : 0.)*(h>epsilon ? Uy/h : 0.) );
}








