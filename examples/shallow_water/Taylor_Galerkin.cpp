#include "Taylor_Galerkin.h"
#include <algorithm>


TG2_scheme::TG2_scheme(const Q1& sol,
                       const Q1& sold,
                       const Q1& soldd,
                       Q1& incr,
                       Q0& sol_onehalf_incr,
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
: sol(sol), sold(sold), soldd(soldd), incr(incr), sol_onehalf_incr(sol_onehalf_incr), 
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

  }
  
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
  
  sol_onehalf_incr[ordh  (quadrant->get_forest_quad_idx ())] = - dt/2.*div_Fh_cell /area;
  sol_onehalf_incr[ordUx (quadrant->get_forest_quad_idx ())] = - dt/2.*div_FUx_cell/area + dt/2.*source_Ux_cell_average;
  sol_onehalf_incr[ordUy (quadrant->get_forest_quad_idx ())] = - dt/2.*div_FUy_cell/area + dt/2.*source_Uy_cell_average;
  
}



void
TG2_scheme::second_step (tmesh::quadrant_iterator quadrant)
{
  
  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;
  
  double h_cell_average = 0., Ux_cell_average = 0., Uy_cell_average = 0.;
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

    h_cell_average  += hdof_c;
    Ux_cell_average += Uxdof_c;
    Uy_cell_average += Uydof_c;

  }
  h_cell_average  /= 4.;
  Ux_cell_average /= 4.;
  Uy_cell_average /= 4.;


  
  // weights coefficients for the flux term
  der_coeffs_x = {-Dy/2.*isdof_or_hanging[0], +Dy/2.*isdof_or_hanging[1],
                  -Dy/2.*isdof_or_hanging[2], +Dy/2.*isdof_or_hanging[3]};
  
  der_coeffs_y = {-Dx/2.*isdof_or_hanging[0], -Dx/2.*isdof_or_hanging[1],
                  +Dx/2.*isdof_or_hanging[2], +Dx/2.*isdof_or_hanging[3]};
  
  double F_star_h_x  = 0., F_star_h_y  = 0.,
         F_star_Ux_x = 0., F_star_Ux_y = 0.,
         F_star_Uy_x = 0., F_star_Uy_y = 0.;
  
  
  const auto & h_der_onehalf  = sol_onehalf_incr[ordh  (quadrant->get_forest_quad_idx ())]; 
  const auto & Ux_der_onehalf = sol_onehalf_incr[ordUx (quadrant->get_forest_quad_idx ())]; 
  const auto & Uy_der_onehalf = sol_onehalf_incr[ordUy (quadrant->get_forest_quad_idx ())]; 
  
  
  // compute local extrema
  // const auto h_min  = *std::min_element(hdof.begin(),  hdof.end() );
  // const auto h_max  = *std::max_element(hdof.begin(),  hdof.end() );
 
  // const auto Ux_min = *std::min_element(Uxdof.begin(), Uxdof.end());
  // const auto Ux_max = *std::max_element(Uxdof.begin(), Uxdof.end());

  // const auto Uy_min = *std::min_element(Uydof.begin(), Uydof.end());
  // const auto Uy_max = *std::max_element(Uydof.begin(), Uydof.end());  

  double h_min  = h_cell_average,  h_max  = h_cell_average;
  double Ux_min = Ux_cell_average, Ux_max = Ux_cell_average;
  double Uy_min = Uy_cell_average, Uy_max = Uy_cell_average;
  for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
       quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
  {

    double h_current_cell = 0., Ux_current_cell = 0., Uy_current_cell = 0.;
    for (int ii = 0; ii < 4; ++ii){
      if (! quadrant_nei->is_hanging (ii)){
        h_current_cell  += sol [ordh  (quadrant_nei->gt (ii))];
        Ux_current_cell += sol [ordUx (quadrant_nei->gt (ii))];
        Uy_current_cell += sol [ordUy (quadrant_nei->gt (ii))];
      } else {
        h_current_cell  += .5 * (sol [ordh  (quadrant_nei->gparent(0,ii))] +
                                 sol [ordh  (quadrant_nei->gparent(1,ii))]);
        Ux_current_cell += .5 * (sol [ordUx (quadrant_nei->gparent(0,ii))] +
                                 sol [ordUx (quadrant_nei->gparent(1,ii))]);
        Uy_current_cell += .5 * (sol [ordUy (quadrant_nei->gparent(0,ii))] +
                                 sol [ordUy (quadrant_nei->gparent(1,ii))]);
      }
    }
    h_current_cell  /= 4.;
    Ux_current_cell /= 4.;
    Uy_current_cell /= 4.;

    h_min  = std::min(h_current_cell,  h_min);
    h_max  = std::max(h_current_cell,  h_max);

    Ux_min = std::min(Ux_current_cell, Ux_min);
    Ux_max = std::max(Ux_current_cell, Ux_max);

    Uy_min = std::min(Uy_current_cell, Uy_min);
    Uy_max = std::max(Uy_current_cell, Uy_max);

  }
  
  for (int ii = 0; ii < 4; ++ii){
    const auto & hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    vel_rusanov_cell_x += hpoint>epsilon ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
    vel_rusanov_cell_y += hpoint>epsilon ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.;
  
  const auto h_cell  = h_cell_average  + h_der_onehalf;
  const auto Ux_cell = Ux_cell_average + Ux_der_onehalf;
  const auto Uy_cell = Uy_cell_average + Uy_der_onehalf;
  
  grad_cell_h  = {.5 * ( (hdof [3] - hdof [2]) + (hdof [1] - hdof [0]) ), .5 * ( (hdof [2] - hdof [0]) + (hdof [3] - hdof [1]) )};
  grad_cell_Ux = {.5 * ( (Uxdof[3] - Uxdof[2]) + (Uxdof[1] - Uxdof[0]) ), .5 * ( (Uxdof[2] - Uxdof[0]) + (Uxdof[3] - Uxdof[1]) )};
  grad_cell_Uy = {.5 * ( (Uydof[3] - Uydof[2]) + (Uydof[1] - Uydof[0]) ), .5 * ( (Uydof[2] - Uydof[0]) + (Uydof[3] - Uydof[1]) )};
  
  // compute flux limiter, grad limiter
  const auto toll = 1e-17;
  double phi_cell_h = 1., phi_cell_Ux = 1., phi_cell_Uy = 1.;
  for (int ii = 0; ii < 4; ++ii){
    
    const auto x_node = quadrant->p(0, ii);
    const auto y_node = quadrant->p(1, ii);
    const auto x_center = quadrant->centroid(0);
    const auto y_center = quadrant->centroid(1);

    const auto Dx_v = x_node - x_center;
    const auto Dy_v = y_node - y_center;

    const auto h_vertex  = h_cell_average  + grad_cell_h [0]*Dx_v/Dx + grad_cell_h [1]*Dy_v/Dy + h_der_onehalf *2.;
    const auto Ux_vertex = Ux_cell_average + grad_cell_Ux[0]*Dx_v/Dx + grad_cell_Ux[1]*Dy_v/Dy + Ux_der_onehalf*2.;
    const auto Uy_vertex = Uy_cell_average + grad_cell_Uy[0]*Dx_v/Dx + grad_cell_Uy[1]*Dy_v/Dy + Uy_der_onehalf*2.;
    // const auto h_vertex  = sol[ordh  (quadrant->gt (ii))];
    // const auto Ux_vertex = sol[ordUx (quadrant->gt (ii))];
    // const auto Uy_vertex = sol[ordUy (quadrant->gt (ii))];

    flux_limiter(h_min,  h_max,  h_vertex,  h_cell_average,  toll, phi_cell_h );
    flux_limiter(Ux_min, Ux_max, Ux_vertex, Ux_cell_average, toll, phi_cell_Ux);
    flux_limiter(Uy_min, Uy_max, Uy_vertex, Uy_cell_average, toll, phi_cell_Uy);
    
  }
  //std::cout << phi_cell_h << " " << phi_cell_Ux << " " << phi_cell_Uy << std::endl;
  // phi_cell_h  = 0.;
  // phi_cell_Ux = 0.;
  // phi_cell_Uy = 0.;

  phi_cell_h  = 1. - phi_cell_h;
  phi_cell_Ux = 1. - phi_cell_Ux;
  phi_cell_Uy = 1. - phi_cell_Uy;


  const auto correction_term_h_x  = grad_cell_h [0] * vel_rusanov_cell_y * phi_cell_h;
  const auto correction_term_h_y  = grad_cell_h [1] * vel_rusanov_cell_x * phi_cell_h;

  const auto correction_term_Ux_x = grad_cell_Ux[0] * vel_rusanov_cell_y * phi_cell_Ux;
  const auto correction_term_Ux_y = grad_cell_Ux[1] * vel_rusanov_cell_x * phi_cell_Ux;

  const auto correction_term_Uy_x = grad_cell_Uy[0] * vel_rusanov_cell_y * phi_cell_Uy;
  const auto correction_term_Uy_y = grad_cell_Uy[1] * vel_rusanov_cell_x * phi_cell_Uy;
  

  F_star_h_x  = h_flux_formula_x (h_cell, Ux_cell, Uy_cell) - correction_term_h_x;
  F_star_h_y  = h_flux_formula_y (h_cell, Ux_cell, Uy_cell) - correction_term_h_y;

  F_star_Ux_x = Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell) - correction_term_Ux_x;
  F_star_Ux_y = Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell) - correction_term_Ux_y;

  F_star_Uy_x = Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell) - correction_term_Uy_x;
  F_star_Uy_y = Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell) - correction_term_Uy_y;

  
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
    const auto Ux_ = der_coeffs_x[ii] * F_star_Ux_x + der_coeffs_y[ii] * F_star_Ux_y + .25*area*isdof_or_hanging[ii] * Ux_src_formula(h_cell, Ux_cell, Uy_cell, partial_x_Z_average);
    const auto Uy_ = der_coeffs_x[ii] * F_star_Uy_x + der_coeffs_y[ii] * F_star_Uy_y + .25*area*isdof_or_hanging[ii] * Uy_src_formula(h_cell, Ux_cell, Uy_cell, partial_y_Z_average);
    

    if (! quadrant->is_hanging (ii)){
      incr [ordh  (quadrant->gt (ii))] += h_;
      incr [ordUx (quadrant->gt (ii))] += Ux_;
      incr [ordUy (quadrant->gt (ii))] += Uy_;
      
      
      auto boundary_idxx = quadrant->ex (ii);
      auto boundary_idxy = quadrant->ey (ii);
      
      
      double F_star_h_x_b  = 0., F_star_h_y_b  = 0.,
             F_star_Ux_x_b = 0., F_star_Ux_y_b = 0., 
             F_star_Uy_x_b = 0., F_star_Uy_y_b = 0.;


      if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY) {
        
        F_star_h_x_b  = h_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
        F_star_h_y_b  = h_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);

        F_star_Ux_x_b = Ux_flux_formula_x(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
        F_star_Ux_y_b = Ux_flux_formula_y(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);

        F_star_Uy_x_b = Uy_flux_formula_x(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
        F_star_Uy_y_b = Uy_flux_formula_y(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
    
        
        auto h_b  = -der_coeffs_x[ii] * F_star_h_x_b  + der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = -der_coeffs_x[ii] * F_star_Ux_x_b + der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = -der_coeffs_x[ii] * F_star_Uy_x_b + der_coeffs_y[ii] * F_star_Uy_y_b;
        
        
        incr [ordh    (quadrant->gt (ii))] += h_b;
        incr [ordUx   (quadrant->gt (ii))] += Ux_b;
        incr [ordUy   (quadrant->gt (ii))] += Uy_b;
        
        
        for (int jj = 0; jj < 4; ++jj){
          if ( quadrant->is_hanging (jj) && (der_coeffs_y[jj]*der_coeffs_y[ii] > 0.) )
          {
            F_star_h_y_b  = h_flux_formula_y  (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
            F_star_Ux_y_b = Ux_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
            F_star_Uy_y_b = Uy_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
            
            h_b  = 2. * der_coeffs_y[jj] * F_star_h_y_b;
            Ux_b = 2. * der_coeffs_y[jj] * F_star_Ux_y_b; 
            Uy_b = 2. * der_coeffs_y[jj] * F_star_Uy_y_b;
            
            incr [ordh  (quadrant->gt (ii))] += h_b;
            incr [ordUx (quadrant->gt (ii))] += Ux_b; 
            incr [ordUy (quadrant->gt (ii))] += Uy_b;
          }
          
        }
        
        
        
        
      }
      
      if (boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY) {

        F_star_h_x_b = h_flux_formula_x  (h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_h_y_b = h_flux_formula_y  (h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);

        F_star_Ux_x_b = Ux_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_Ux_y_b = Ux_flux_formula_y(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
 
        F_star_Uy_x_b = Uy_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_Uy_y_b = Uy_flux_formula_y(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        
        
        auto h_b  = der_coeffs_x[ii] * F_star_h_x_b  - der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = der_coeffs_x[ii] * F_star_Ux_x_b - der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = der_coeffs_x[ii] * F_star_Uy_x_b - der_coeffs_y[ii] * F_star_Uy_y_b;
        

        incr [ordh    (quadrant->gt (ii))] += h_b;
        incr [ordUx   (quadrant->gt (ii))] += Ux_b;
        incr [ordUy   (quadrant->gt (ii))] += Uy_b;
        
        
        for (int jj = 0; jj < 4; ++jj){
          if ( quadrant->is_hanging (jj) && (der_coeffs_x[jj]*der_coeffs_x[ii] > 0.) )
          {
            F_star_h_x_b  = h_flux_formula_x (h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
            F_star_Ux_x_b = Ux_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
            F_star_Uy_x_b = Uy_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
            
            h_b  = 2. * der_coeffs_x[jj] * F_star_h_x_b;
            Ux_b = 2. * der_coeffs_x[jj] * F_star_Ux_x_b;
            Uy_b = 2. * der_coeffs_x[jj] * F_star_Uy_x_b;
            
            
            incr [ordh    (quadrant->gt (ii))] += h_b;
            incr [ordUx   (quadrant->gt (ii))] += Ux_b;
            incr [ordUy   (quadrant->gt (ii))] += Uy_b;
          }
          
        }
        
        
      }
      
      // corner points!
      if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY && boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY)
      {
        
        F_star_h_x_b  = h_flux_formula_x  (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_h_y_b  = h_flux_formula_y  (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        
        F_star_Ux_x_b = Ux_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_Ux_y_b = Ux_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        
        F_star_Uy_x_b = Uy_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_Uy_y_b = Uy_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        
        
        auto h_b  = -der_coeffs_x[ii] * F_star_h_x_b  - der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = -der_coeffs_x[ii] * F_star_Ux_x_b - der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = -der_coeffs_x[ii] * F_star_Uy_x_b - der_coeffs_y[ii] * F_star_Uy_y_b;
        

        incr [ordh    (quadrant->gt (ii))] += h_b;
        incr [ordUx   (quadrant->gt (ii))] += Ux_b;
        incr [ordUy   (quadrant->gt (ii))] += Uy_b;
        
      }
      
      
    } else {
      
      incr [ordh  (quadrant->gparent(0,ii))] += h_;
      incr [ordh  (quadrant->gparent(1,ii))] += h_;
      
      incr [ordUx (quadrant->gparent(0,ii))] += Ux_;
      incr [ordUx (quadrant->gparent(1,ii))] += Ux_;
      
      incr [ordUy (quadrant->gparent(0,ii))] += Uy_;
      incr [ordUy (quadrant->gparent(1,ii))] += Uy_;
      
    }
  }
}

void
TG2_scheme::flux_limiter(const double& Q_min, const double& Q_max, const double& Q_vertex, const double& Q_cell, const double& toll, double& phi_cell_Q)
{
  
  if (std::abs(Q_vertex-Q_cell) <= toll)
  {
    phi_cell_Q = std::min(1., phi_cell_Q);
  }
  else if ((Q_vertex-Q_cell) < -toll)
  {
    phi_cell_Q = std::min(std::min(1., (Q_min-Q_cell)/(Q_vertex-Q_cell)), phi_cell_Q );
  }
  else //((Q_vertex-Q_cell) > toll)
  {
    phi_cell_Q = std::min(std::min(1., (Q_max-Q_cell)/(Q_vertex-Q_cell)), phi_cell_Q );
  }
  
  
  if (Q_cell == 0 || Q_vertex == 0)
  {
    phi_cell_Q = 0.; 
  }
  
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
{ return Ux; }

double
TG2_scheme::h_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return Uy; }

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



