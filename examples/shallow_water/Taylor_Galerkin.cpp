#include "Taylor_Galerkin.h"
#include <algorithm>
#include <cassert>

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
                       const Q0& slope_x,
                       const Q0& slope_y,
                       const double& DELTAT,
                       const double& h_min,
                       const bool& is_non_reflBC,
                       const bool& is_bed_friction,
                       const bool& is_stress_tensor,
                       const double& grav,
                       const double& density,
                       const double& turbulence_coeff,
                       const double& surface_pressure, 
                       const double& bed_friction_angle_rad,
                       const double& fluid_viscosity,
                       const double& yield_shear_stress)
: sol(sol), sold(sold), soldd(soldd), incr(incr), incr_anti_diff(incr_anti_diff), P_plus(P_plus), P_minus(P_minus), sol_onehalf(sol_onehalf), mass(mass),
  ordh(oh), ordUx(oUx), ordUy(oUy), Z(Z), slope_x(slope_x), slope_y(slope_y), DELTAT(DELTAT), epsilon(h_min), is_non_reflBC(is_non_reflBC), is_bed_friction(is_bed_friction), is_stress_tensor(is_stress_tensor), grav(grav),
  density(density), turbulence_coeff(turbulence_coeff), surface_pressure(surface_pressure), bed_friction_angle_rad(bed_friction_angle_rad), fluid_viscosity(fluid_viscosity), yield_shear_stress(yield_shear_stress)
{ }

 

void
TG2_scheme::compute_dt (tmesh::quadrant_iterator quadrant)
{

  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 

  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii); 
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];

  std::array<double,4> vel_x, vel_y;  
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

    const auto & hpoint = hdof[ii];
    vel_x[ii] = hpoint>epsilon ? Uxdof[ii]/hpoint : 0.;
    vel_y[ii] = hpoint>epsilon ? Uydof[ii]/hpoint : 0.;
  }
  
  grad_cell_ux = {.5 * ( (vel_x[3] - vel_x[2]) + (vel_x[1] - vel_x[0]) )/Dx, .5 * ( (vel_x[2] - vel_x[0]) + (vel_x[3] - vel_x[1]) )/Dy};
  grad_cell_uy = {.5 * ( (vel_y[3] - vel_y[2]) + (vel_y[1] - vel_y[0]) )/Dx, .5 * ( (vel_y[2] - vel_y[0]) + (vel_y[3] - vel_y[1]) )/Dy};
  
  for (int ii = 0; ii < 4; ++ii){
    const auto & hpoint  = hdof [ii];
    const auto & Uxpoint = Uxdof[ii];
    const auto & Uypoint = Uydof[ii];

    const auto & vel_x_point = vel_x[ii];
    const auto & vel_y_point = vel_y[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    
    const auto vel_rusanov_x = hpoint>epsilon ? std::abs(vel_x_point)+celerity : 0.;
    const auto vel_rusanov_y = hpoint>epsilon ? std::abs(vel_y_point)+celerity : 0.;
    
    const auto dtx_adv = hpoint>epsilon ? Dx/vel_rusanov_x : dt; //dt; // DELTAT;
    const auto dty_adv = hpoint>epsilon ? Dy/vel_rusanov_y : dt; //dt; // DELTAT;

    // std::array<double,6> def_grad = compute_nodal_def_grad (hpoint, Uxpoint, Uypoint, grad_cell_ux, grad_cell_uy);

    // double second_invariant = 0.;
    // for (int i_def = 0; i_def < 6; i_def++)
    // {
    //   second_invariant += def_grad[i_def]*def_grad[i_def];
    // }
    // second_invariant *= .5;

    // const double viscos = second_invariant!=0 ? .5*(yield_shear_stress/std::sqrt(second_invariant) + 2*fluid_viscosity)/density : 0.;

    // const auto dtx_dif = viscos!=0 ? Dx*Dx*.5/viscos : dt; //dt; // DELTAT;
    // const auto dty_dif = viscos!=0 ? Dy*Dy*.5/viscos : dt; //dt; // DELTAT;

    const auto dtoptx = dtx_adv; //std::min( dtx_adv, 2.*std::abs(Uxpoint)/std::abs(Ux_src_formula(hpoint, Uxpoint, Uypoint, slope_x[index_quadrant])) );///(dtx_adv + dtx_dif); //hpoint>epsilon ? Dx/vel_rusanov_x : dt;
    const auto dtopty = dty_adv; //std::min( dty_adv, 2.*std::abs(Uypoint)/std::abs(Uy_src_formula(hpoint, Uxpoint, Uypoint, slope_y[index_quadrant])) );//(dty_adv * dty_dif)/(dty_adv + dty_dif);///(dty_adv + dty_dif); //hpoint>epsilon ? Dy/vel_rusanov_y : dt;
    const auto dtopt = dtoptx > dtopty ? dtopty : dtoptx;
    if (dt > dtopt) set_dt (dtopt);
    
    //std::cout << dt << std::endl;
    //if (dtx_dif!=dt || dty_dif!=dt)

    // if (dt<=0.)
    // {
    //   std::cout << dt << " " << viscos << " " << dtx_dif << " " << dty_dif << " " << dtx_adv << " " << dty_adv << " " << 2.*Uypoint/std::abs(Uy_src_formula(hpoint, Uxpoint, Uypoint, slope_y[index_quadrant])) << std::endl;
    //   exit(1);
    // }
    

    
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
  nu_htot += Nu_hmean_cell*(time-timed);
  
}

void
TG2_scheme::first_step (tmesh::quadrant_iterator quadrant)
{
  
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 
  
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii);
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;


  
  double h_cell_average = 0., Ux_cell_average = 0., Uy_cell_average = 0., source_Ux_cell_average = 0., source_Uy_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    double hdof_c, Uxdof_c, Uydof_c;

    if (! quadrant->is_hanging (ii) )
    {
      hdof_c    = sol [ordh    (quadrant->gt (ii) )];
      Uxdof_c   = sol [ordUx   (quadrant->gt (ii) )];
      Uydof_c   = sol [ordUy   (quadrant->gt (ii) )];
    }
    else
    {
      hdof_c    = .5 * (sol [ordh    (quadrant->gparent (0, ii) )] +
                        sol [ordh    (quadrant->gparent (1, ii) )]);
      Uxdof_c   = .5 * (sol [ordUx   (quadrant->gparent (0, ii) )] +
                        sol [ordUx   (quadrant->gparent (1, ii) )]);
      Uydof_c   = .5 * (sol [ordUy   (quadrant->gparent (0, ii) )] +
                        sol [ordUy   (quadrant->gparent (1, ii) )]);
    }

    h_cell_average  += hdof_c;
    Ux_cell_average += Uxdof_c;
    Uy_cell_average += Uydof_c;

    source_Ux_cell_average += Ux_src_formula (hdof_c, Uxdof_c, Uydof_c, slope_x[index_quadrant]); 
    source_Uy_cell_average += Uy_src_formula (hdof_c, Uxdof_c, Uydof_c, slope_y[index_quadrant]);
    
    fluxx_h_node[ii]  = h_flux_formula_x   (hdof_c, Uxdof_c, Uydof_c);
    fluxy_h_node[ii]  = h_flux_formula_y   (hdof_c, Uxdof_c, Uydof_c);
    fluxx_Ux_node[ii] = Ux_flux_formula_x  (hdof_c, Uxdof_c, Uydof_c);
    fluxy_Ux_node[ii] = Ux_flux_formula_y  (hdof_c, Uxdof_c, Uydof_c);
    fluxx_Uy_node[ii] = Uy_flux_formula_x  (hdof_c, Uxdof_c, Uydof_c);
    fluxy_Uy_node[ii] = Uy_flux_formula_y  (hdof_c, Uxdof_c, Uydof_c);

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
  
  sol_onehalf[ordh  (index_quadrant)] = h_cell_average  - dt / 2. * div_Fh_cell  / area;
  sol_onehalf[ordUx (index_quadrant)] = Ux_cell_average - dt / 2. * div_FUx_cell / area + dt / 2. * source_Ux_cell_average;
  sol_onehalf[ordUy (index_quadrant)] = Uy_cell_average - dt / 2. * div_FUy_cell / area + dt / 2. * source_Uy_cell_average;
  


}


void
TG2_scheme::compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant)
{

  // look at tmesh.h
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 
  
  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;

  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){
      hdof[ii]    = sol [ordh    (quadrant->gt (ii))];
      Uxdof[ii]   = sol [ordUx   (quadrant->gt (ii))];
      Uydof[ii]   = sol [ordUy   (quadrant->gt (ii))];
      
      isdof_or_hanging[ii] = 1.;
    } else {
      hdof[ii]    = .5 * (sol [ordh  (quadrant->gparent(0,ii))] +
                          sol [ordh  (quadrant->gparent(1,ii))]);
      Uxdof[ii]   = .5 * (sol [ordUx (quadrant->gparent(0,ii))] +
                          sol [ordUx (quadrant->gparent(1,ii))]);
      Uydof[ii]   = .5 * (sol [ordUy (quadrant->gparent(0,ii))] +
                          sol [ordUy (quadrant->gparent(1,ii))]);
      
      isdof_or_hanging[ii] = .5;
    }
  }

  // weights coefficients for the flux term
  der_coeffs_x = {-Dy/2.*isdof_or_hanging[0], +Dy/2.*isdof_or_hanging[1],
    -Dy/2.*isdof_or_hanging[2], +Dy/2.*isdof_or_hanging[3]};
  
  der_coeffs_y = {-Dx/2.*isdof_or_hanging[0], -Dx/2.*isdof_or_hanging[1],
    +Dx/2.*isdof_or_hanging[2], +Dx/2.*isdof_or_hanging[3]};

  std::array<double,4> vel_x, vel_y;
  double vel_rusanov_cell_x = 0., vel_rusanov_cell_y = 0.;
  for (int ii = 0; ii < 4; ++ii){
    const auto& hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    vel_x[ii] = hpoint>epsilon ? Uxdof[ii]/hpoint : 0.;
    vel_y[ii] = hpoint>epsilon ? Uydof[ii]/hpoint : 0.;
    vel_rusanov_cell_x += hpoint>epsilon ? std::abs(vel_x[ii])+celerity : 0.;
    vel_rusanov_cell_y += hpoint>epsilon ? std::abs(vel_y[ii])+celerity : 0.;
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.; 

  grad_cell_h  = {.5 * ( (hdof [3] - hdof [2]) + (hdof [1] - hdof [0]) )   , .5 * ( (hdof [2] - hdof [0]) + (hdof [3] - hdof [1]) )};
  grad_cell_Ux = {.5 * ( (Uxdof[3] - Uxdof[2]) + (Uxdof[1] - Uxdof[0]) )   , .5 * ( (Uxdof[2] - Uxdof[0]) + (Uxdof[3] - Uxdof[1]) )};
  grad_cell_Uy = {.5 * ( (Uydof[3] - Uydof[2]) + (Uydof[1] - Uydof[0]) )   , .5 * ( (Uydof[2] - Uydof[0]) + (Uydof[3] - Uydof[1]) )};

  grad_cell_ux = {.5 * ( (vel_x[3] - vel_x[2]) + (vel_x[1] - vel_x[0]) )/Dx, .5 * ( (vel_x[2] - vel_x[0]) + (vel_x[3] - vel_x[1]) )/Dy};
  grad_cell_uy = {.5 * ( (vel_y[3] - vel_y[2]) + (vel_y[1] - vel_y[0]) )/Dx, .5 * ( (vel_y[2] - vel_y[0]) + (vel_y[3] - vel_y[1]) )/Dy};
  

  const double & h_cell  = sol_onehalf[ordh  (index_quadrant)];
  const double & Ux_cell = sol_onehalf[ordUx (index_quadrant)];
  const double & Uy_cell = sol_onehalf[ordUy (index_quadrant)];


  const auto diff_term_h_x  = grad_cell_h [0] * vel_rusanov_cell_y;
  const auto diff_term_h_y  = grad_cell_h [1] * vel_rusanov_cell_x;

  const auto diff_term_Ux_x = grad_cell_Ux[0] * vel_rusanov_cell_y;
  const auto diff_term_Ux_y = grad_cell_Ux[1] * vel_rusanov_cell_x;

  const auto diff_term_Uy_x = grad_cell_Uy[0] * vel_rusanov_cell_y;
  const auto diff_term_Uy_y = grad_cell_Uy[1] * vel_rusanov_cell_x;
  


  const auto F_star_h_x  = h_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_h_x;
  const auto F_star_h_y  = h_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_h_y;

  const auto F_star_Ux_x = Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_x;
  const auto F_star_Ux_y = Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_y;

  const auto F_star_Uy_x = Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_x;
  const auto F_star_Uy_y = Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_y;






  for (int ii = 0; ii < 4; ++ii){
    

    sigma_stress = is_stress_tensor ? compute_nodal_stress (hdof[ii], Uxdof[ii], Uydof[ii], grad_cell_ux, grad_cell_uy) : std::array<double,3>{{0., 0., 0.}};
    

    const double D_Ux_x = is_stress_tensor ? Ux_stress_formula_x(hdof[ii], Uxdof[ii], Uydof[ii]) : 0.;
    const double D_Ux_y = is_stress_tensor ? Ux_stress_formula_y(hdof[ii], Uxdof[ii], Uydof[ii]) : 0.;

    const double D_Uy_x = is_stress_tensor ? Uy_stress_formula_x(hdof[ii], Uxdof[ii], Uydof[ii]) : 0.;
    const double D_Uy_y = is_stress_tensor ? Uy_stress_formula_y(hdof[ii], Uxdof[ii], Uydof[ii]) : 0.;

    if (std::isnan(D_Ux_x) || std::isnan(D_Ux_y) || std::isnan(D_Uy_x) || std::isnan(D_Uy_y))
    {
      std::cout << D_Ux_x << " " << D_Ux_y << " " << D_Uy_x << " " << D_Uy_y << std::endl;
      exit(1);
    }
    

    const auto h_  = der_coeffs_x[ii]*F_star_h_x +der_coeffs_y[ii]*F_star_h_y;
    const auto Ux_ = der_coeffs_x[ii]*F_star_Ux_x+der_coeffs_y[ii]*F_star_Ux_y + der_coeffs_x[ii]*(1./3.)*D_Ux_x+der_coeffs_y[ii]*(1./3.)*D_Ux_y + .25*area*isdof_or_hanging[ii]*Ux_src_formula(h_cell, Ux_cell, Uy_cell, slope_x[index_quadrant]);
    const auto Uy_ = der_coeffs_x[ii]*F_star_Uy_x+der_coeffs_y[ii]*F_star_Uy_y + der_coeffs_x[ii]*(1./3.)*D_Uy_x+der_coeffs_y[ii]*(1./3.)*D_Uy_y + .25*area*isdof_or_hanging[ii]*Uy_src_formula(h_cell, Ux_cell, Uy_cell, slope_y[index_quadrant]);
    
    const auto h_al  = der_coeffs_x[ii] * diff_term_h_x  + der_coeffs_y[ii] * diff_term_h_y;
    const auto Ux_al = der_coeffs_x[ii] * diff_term_Ux_x + der_coeffs_y[ii] * diff_term_Ux_y;
    const auto Uy_al = der_coeffs_x[ii] * diff_term_Uy_x + der_coeffs_y[ii] * diff_term_Uy_y;

    incr_anti_diff[ordh (index_quadrant)][ii] = h_al;
    incr_anti_diff[ordUx(index_quadrant)][ii] = Ux_al;
    incr_anti_diff[ordUy(index_quadrant)][ii] = Uy_al;


    if (! quadrant->is_hanging (ii)){

      incr [ordh  (quadrant->gt (ii))] += h_;
      incr [ordUx (quadrant->gt (ii))] += Ux_;
      incr [ordUy (quadrant->gt (ii))] += Uy_;

      P_plus [ordh  (quadrant->gt (ii))] += std::max(0., h_al );
      P_plus [ordUx (quadrant->gt (ii))] += std::max(0., Ux_al);
      P_plus [ordUy (quadrant->gt (ii))] += std::max(0., Uy_al);

      P_minus [ordh  (quadrant->gt (ii))] += std::min(0., h_al );
      P_minus [ordUx (quadrant->gt (ii))] += std::min(0., Ux_al);
      P_minus [ordUy (quadrant->gt (ii))] += std::min(0., Uy_al);
      
      
      const auto boundary_idxx = quadrant->ex (ii);
      const auto boundary_idxy = quadrant->ey (ii);
      
      
      
      if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY) {
        double F_star_h_x_b = 0., F_star_h_y_b = 0.,
        F_star_Ux_x_b = 0., F_star_Ux_y_b = 0., F_star_Uy_x_b = 0., F_star_Uy_y_b = 0.;
        
        F_star_h_x_b  = h_flux_formula_x(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
        F_star_h_y_b  = h_flux_formula_y(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);

        F_star_Ux_x_b = Ux_flux_formula_x(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
        F_star_Ux_y_b = Ux_flux_formula_y(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);

        F_star_Uy_x_b = Uy_flux_formula_x(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
        F_star_Uy_y_b = Uy_flux_formula_y(h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, Uy_cell);
    
        
        
        auto h_b  = -der_coeffs_x[ii]*F_star_h_x_b  + der_coeffs_y[ii]*F_star_h_y_b;
        auto Ux_b = -der_coeffs_x[ii]*F_star_Ux_x_b + der_coeffs_y[ii]*F_star_Ux_y_b;
        auto Uy_b = -der_coeffs_x[ii]*F_star_Uy_x_b + der_coeffs_y[ii]*F_star_Uy_y_b;

        auto h_a  = -der_coeffs_x[ii]*diff_term_h_x  + der_coeffs_y[ii]*diff_term_h_y;
        auto Ux_a = -der_coeffs_x[ii]*diff_term_Ux_x + der_coeffs_y[ii]*diff_term_Ux_y;
        auto Uy_a = -der_coeffs_x[ii]*diff_term_Uy_x + der_coeffs_y[ii]*diff_term_Uy_y;
        
        
        
        incr [ordh  (quadrant->gt (ii))] += h_b;
        incr [ordUx (quadrant->gt (ii))] += Ux_b;
        incr [ordUy (quadrant->gt (ii))] += Uy_b;

        incr_anti_diff[ordh  (index_quadrant)][ii] += h_a;
        incr_anti_diff[ordUx (index_quadrant)][ii] += Ux_a;
        incr_anti_diff[ordUy (index_quadrant)][ii] += Uy_a;

        P_plus [ordh  (quadrant->gt (ii))] += std::max(0., h_a );
        P_plus [ordUx (quadrant->gt (ii))] += std::max(0., Ux_a);
        P_plus [ordUy (quadrant->gt (ii))] += std::max(0., Uy_a);

        P_minus [ordh  (quadrant->gt (ii))] += std::min(0., h_a );
        P_minus [ordUx (quadrant->gt (ii))] += std::min(0., Ux_a);
        P_minus [ordUy (quadrant->gt (ii))] += std::min(0., Uy_a);
        
        
        for (int jj = 0; jj < 4; ++jj){
          if ( quadrant->is_hanging (jj) && (der_coeffs_y[jj]*der_coeffs_y[ii] > 0.) )
          {
            
            h_b  = 2. * der_coeffs_y[jj] * F_star_h_y_b;
            Ux_b = 2. * der_coeffs_y[jj] * F_star_Ux_y_b;
            Uy_b = 2. * der_coeffs_y[jj] * F_star_Uy_y_b;

            h_a  = 2. * der_coeffs_y[jj] * diff_term_h_y;
            Ux_a = 2. * der_coeffs_y[jj] * diff_term_Ux_y;
            Uy_a = 2. * der_coeffs_y[jj] * diff_term_Uy_y;
            
            
            incr [ordh  (quadrant->gt (ii))] += h_b;
            incr [ordUx (quadrant->gt (ii))] += Ux_b;
            incr [ordUy (quadrant->gt (ii))] += Uy_b;

            incr_anti_diff[ordh  (index_quadrant)][ii] += h_a;
            incr_anti_diff[ordUx (index_quadrant)][ii] += Ux_a;
            incr_anti_diff[ordUy (index_quadrant)][ii] += Uy_a;
            
            P_plus [ordh  (quadrant->gt (ii))] += std::max(0., h_a );
            P_plus [ordUx (quadrant->gt (ii))] += std::max(0., Ux_a); 
            P_plus [ordUy (quadrant->gt (ii))] += std::max(0., Uy_a);

            P_minus [ordh  (quadrant->gt (ii))] += std::min(0., h_a );
            P_minus [ordUx (quadrant->gt (ii))] += std::min(0., Ux_a);
            P_minus [ordUy (quadrant->gt (ii))] += std::min(0., Uy_a);
          }
          
        }
        
        
        
        
      }
      
      if (boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY) {
        double F_star_h_x_b = 0., F_star_h_y_b = 0.,
        F_star_Ux_x_b = 0., F_star_Ux_y_b = 0., F_star_Uy_x_b = 0., F_star_Uy_y_b = 0.;
        
        F_star_h_x_b  = h_flux_formula_x (h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_h_y_b  = h_flux_formula_y (h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);

        F_star_Ux_x_b = Ux_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_Ux_y_b = Ux_flux_formula_y(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);

        F_star_Uy_x_b = Uy_flux_formula_x(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_Uy_y_b = Uy_flux_formula_y(h_cell, Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        
        
        auto h_b    = der_coeffs_x[ii] * F_star_h_x_b    - der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b   = der_coeffs_x[ii] * F_star_Ux_x_b   - der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b   = der_coeffs_x[ii] * F_star_Uy_x_b   - der_coeffs_y[ii] * F_star_Uy_y_b;

        auto h_a    = der_coeffs_x[ii] * diff_term_h_x   - der_coeffs_y[ii] * diff_term_h_y;
        auto Ux_a   = der_coeffs_x[ii] * diff_term_Ux_x  - der_coeffs_y[ii] * diff_term_Ux_y;
        auto Uy_a   = der_coeffs_x[ii] * diff_term_Uy_x  - der_coeffs_y[ii] * diff_term_Uy_y;
              
        
        incr [ordh    (quadrant->gt (ii))] += h_b;
        incr [ordUx   (quadrant->gt (ii))] += Ux_b;
        incr [ordUy   (quadrant->gt (ii))] += Uy_b;
        
        incr_anti_diff[ordh (index_quadrant)][ii] += h_a;
        incr_anti_diff[ordUx(index_quadrant)][ii] += Ux_a;
        incr_anti_diff[ordUy(index_quadrant)][ii] += Uy_a;

        P_plus [ordh    (quadrant->gt (ii))] += std::max(0., h_a );
        P_plus [ordUx   (quadrant->gt (ii))] += std::max(0., Ux_a);
        P_plus [ordUy   (quadrant->gt (ii))] += std::max(0., Uy_a);

        P_minus [ordh    (quadrant->gt (ii))] += std::min(0., h_a );
        P_minus [ordUx   (quadrant->gt (ii))] += std::min(0., Ux_a);
        P_minus [ordUy   (quadrant->gt (ii))] += std::min(0., Uy_a);


        for (int jj = 0; jj < 4; ++jj){
          if ( quadrant->is_hanging (jj) && (der_coeffs_x[jj]*der_coeffs_x[ii] > 0.) )
          {
            
            h_b    = 2. * der_coeffs_x[jj] * F_star_h_x_b;
            Ux_b   = 2. * der_coeffs_x[jj] * F_star_Ux_x_b;
            Uy_b   = 2. * der_coeffs_x[jj] * F_star_Uy_x_b;

            h_a    = 2. * der_coeffs_x[jj] * diff_term_h_x;
            Ux_a   = 2. * der_coeffs_x[jj] * diff_term_Ux_x;
            Uy_a   = 2. * der_coeffs_x[jj] * diff_term_Uy_x;
            
            
            incr [ordh    (quadrant->gt (ii))] += h_b;
            incr [ordUx   (quadrant->gt (ii))] += Ux_b;
            incr [ordUy   (quadrant->gt (ii))] += Uy_b;

            incr_anti_diff[ordh (index_quadrant)][ii] += h_a;
            incr_anti_diff[ordUx(index_quadrant)][ii] += Ux_a;
            incr_anti_diff[ordUy(index_quadrant)][ii] += Uy_a;

            P_plus [ordh    (quadrant->gt (ii))] += std::max(0., h_a );
            P_plus [ordUx   (quadrant->gt (ii))] += std::max(0., Ux_a);
            P_plus [ordUy   (quadrant->gt (ii))] += std::max(0., Uy_a);

            P_minus [ordh    (quadrant->gt (ii))] += std::min(0., h_a );
            P_minus [ordUx   (quadrant->gt (ii))] += std::min(0., Ux_a);
            P_minus [ordUy   (quadrant->gt (ii))] += std::min(0., Uy_a);
          }
          
        }
        
        
      }
      
      // corner points!
      if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY && boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY)
      {
        double F_star_h_x_b = 0., F_star_h_y_b = 0.,
        F_star_Ux_x_b = 0., F_star_Ux_y_b = 0., F_star_Uy_x_b = 0., F_star_Uy_y_b = 0.;
        
        F_star_h_x_b = h_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_h_y_b = h_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        
        F_star_Ux_x_b = Ux_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_Ux_y_b = Ux_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        
        F_star_Uy_x_b = Uy_flux_formula_x (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        F_star_Uy_y_b = Uy_flux_formula_y (h_cell, is_non_reflBC ? Ux_cell : -Ux_cell, is_non_reflBC ? Uy_cell : -Uy_cell);
        
        
        auto h_b    = -der_coeffs_x[ii] * F_star_h_x_b    - der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b   = -der_coeffs_x[ii] * F_star_Ux_x_b   - der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b   = -der_coeffs_x[ii] * F_star_Uy_x_b   - der_coeffs_y[ii] * F_star_Uy_y_b;

        auto h_a    = -der_coeffs_x[ii] * diff_term_h_x   - der_coeffs_y[ii] * diff_term_h_y;
        auto Ux_a   = -der_coeffs_x[ii] * diff_term_Ux_x  - der_coeffs_y[ii] * diff_term_Ux_y;
        auto Uy_a   = -der_coeffs_x[ii] * diff_term_Uy_x  - der_coeffs_y[ii] * diff_term_Uy_y;
        
        
        
        incr [ordh    (quadrant->gt (ii))] += h_b;
        incr [ordUx   (quadrant->gt (ii))] += Ux_b;
        incr [ordUy   (quadrant->gt (ii))] += Uy_b;

        incr_anti_diff[ordh  (index_quadrant)][ii] += h_a;
        incr_anti_diff[ordUx (index_quadrant)][ii] += Ux_a;
        incr_anti_diff[ordUy (index_quadrant)][ii] += Uy_a;

        P_plus [ordh    (quadrant->gt (ii))] += std::max(0., h_a );
        P_plus [ordUx   (quadrant->gt (ii))] += std::max(0., Ux_a);
        P_plus [ordUy   (quadrant->gt (ii))] += std::max(0., Uy_a);

        P_minus [ordh    (quadrant->gt (ii))] += std::min(0., h_a );
        P_minus [ordUx   (quadrant->gt (ii))] += std::min(0., Ux_a);
        P_minus [ordUy   (quadrant->gt (ii))] += std::min(0., Uy_a);
        
      }
      
      
    } else {
      
      incr [ordh  (quadrant->gparent(0,ii))] += h_;
      incr [ordh  (quadrant->gparent(1,ii))] += h_;
      
      incr [ordUx (quadrant->gparent(0,ii))] += Ux_;
      incr [ordUx (quadrant->gparent(1,ii))] += Ux_;
      
      incr [ordUy (quadrant->gparent(0,ii))] += Uy_;
      incr [ordUy (quadrant->gparent(1,ii))] += Uy_;


      P_plus [ordh  (quadrant->gparent(0,ii))] += std::max(0., h_al);
      P_plus [ordh  (quadrant->gparent(1,ii))] += std::max(0., h_al);
      
      P_plus [ordUx (quadrant->gparent(0,ii))] += std::max(0., Ux_al);
      P_plus [ordUx (quadrant->gparent(1,ii))] += std::max(0., Ux_al);
      
      P_plus [ordUy (quadrant->gparent(0,ii))] += std::max(0., Uy_al);
      P_plus [ordUy (quadrant->gparent(1,ii))] += std::max(0., Uy_al);



      P_minus [ordh  (quadrant->gparent(0,ii))] += std::min(0., h_al);
      P_minus [ordh  (quadrant->gparent(1,ii))] += std::min(0., h_al);
      
      P_minus [ordUx (quadrant->gparent(0,ii))] += std::min(0., Ux_al);
      P_minus [ordUx (quadrant->gparent(1,ii))] += std::min(0., Ux_al);
      
      P_minus [ordUy (quadrant->gparent(0,ii))] += std::min(0., Uy_al);
      P_minus [ordUy (quadrant->gparent(1,ii))] += std::min(0., Uy_al);
      
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



  
  // compute flux correction
  double phi_cell_h = 1., phi_cell_Ux = 1., phi_cell_Uy = 1.;
  for (int ii = 0; ii < 4; ++ii){

    const auto & flux_on_the_node_h  = incr_anti_diff[ordh (index_quadrant)][ii];
    const auto & flux_on_the_node_Ux = incr_anti_diff[ordUx(index_quadrant)][ii];
    const auto & flux_on_the_node_Uy = incr_anti_diff[ordUy(index_quadrant)][ii];


    const auto & hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    const auto vel_rusanov_cell_x = hpoint>epsilon ? (std::abs(Uxdof[ii]/hpoint)+celerity) : 0.;
    const auto vel_rusanov_cell_y = hpoint>epsilon ? (std::abs(Uydof[ii]/hpoint)+celerity) : 0.;

    const auto vel_square_rusanov_cell = vel_rusanov_cell_x * vel_rusanov_cell_y;

    flux_limiter(h_min [ii], h_max [ii], hdof [ii], P_plus_h_dof [ii], P_minus_h_dof  [ii], flux_on_the_node_h,  vel_square_rusanov_cell, phi_cell_h );
    flux_limiter(Ux_min[ii], Ux_max[ii], Uxdof[ii], P_plus_Ux_dof[ii], P_minus_Ux_dof [ii], flux_on_the_node_Ux, vel_square_rusanov_cell, phi_cell_Ux);
    flux_limiter(Uy_min[ii], Uy_max[ii], Uydof[ii], P_plus_Uy_dof[ii], P_minus_Uy_dof [ii], flux_on_the_node_Uy, vel_square_rusanov_cell, phi_cell_Uy);
  }



  for (int ii = 0; ii < 4; ++ii){
 
    const auto flux_on_the_node_h  = incr_anti_diff[ordh (index_quadrant)][ii]*phi_cell_h;
    const auto flux_on_the_node_Ux = incr_anti_diff[ordUx(index_quadrant)][ii]*phi_cell_Ux;
    const auto flux_on_the_node_Uy = incr_anti_diff[ordUy(index_quadrant)][ii]*phi_cell_Uy;

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
TG2_scheme::flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& vel_square_rusanov_cell, double& phi_cell_Q)
{
  const auto Q_plus  = (Q_max-Q_dof)*dt*vel_square_rusanov_cell;
  const auto Q_minus = (Q_min-Q_dof)*dt*vel_square_rusanov_cell;

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
{ return Ux; }

double
TG2_scheme::h_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return Uy; }

double
TG2_scheme::Ux_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ return (h>epsilon ? Ux*Ux/h + grav*h*h/2. : 0.); }
 
double
TG2_scheme::Ux_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return (h>epsilon ? Uy*Ux/h : 0.); }

double
TG2_scheme::Uy_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ return (h>epsilon ? Uy*Ux/h : 0.); }

double
TG2_scheme::Uy_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return (h>epsilon ? Uy*Uy/h + grav*h*h/2. : 0.); }


// stress functions
double
TG2_scheme::Ux_stress_formula_x (const double& h, const double& Ux, const double& Uy)
{ return (-sigma_stress[0]*h/density); }

double
TG2_scheme::Ux_stress_formula_y (const double& h, const double& Ux, const double& Uy)
{ return (-sigma_stress[2]*h/density); }

double
TG2_scheme::Uy_stress_formula_x (const double& h, const double& Ux, const double& Uy)
{ return (-sigma_stress[2]*h/density); }

double
TG2_scheme::Uy_stress_formula_y (const double& h, const double& Ux, const double& Uy)
{ return (-sigma_stress[1]*h/density); }


std::array<double,3>
TG2_scheme::compute_nodal_stress (const double& h, const double& Ux, const double& Uy, const std::array<double,2>& grad_cell_ux, const std::array<double,2>& grad_cell_uy)
{
  // compute \sigma_xx, ...

  // def_grad = [D11, D22, D33, D12, D23, D31]
  // sigma = [sigma_11, sigma_22, sigma_12]

  std::array<double,6> def_grad = compute_nodal_def_grad (h, Ux, Uy, grad_cell_ux, grad_cell_uy);


  double second_invariant = 0.;
  for (int i_def = 0; i_def < 6; i_def++)
  {
    second_invariant += def_grad[i_def]*def_grad[i_def];
  }
  second_invariant *= .5;

  const double viscos = second_invariant!=0 ? yield_shear_stress/std::sqrt(second_invariant) + 2*fluid_viscosity : 0.;

  if (std::isnan(def_grad[0]) || std::isnan(def_grad[1]) || std::isnan(def_grad[2]) || std::isnan(viscos))
  {
    std::cout << grad_cell_ux[0] << " " << grad_cell_ux[1] << " " << viscos << " " << second_invariant << " " << def_grad[0] << " " << def_grad[1] << " " << def_grad[2] << " " << def_grad[3] << " " << def_grad[4] << " " << def_grad[5] << std::endl;
    exit(1);
  }
  //std::cout << viscos << " " << second_invariant << " " << def_grad[0] << " " << def_grad[1] << " " << def_grad[2] << std::endl;
  return(std::array<double,3>{{viscos*def_grad[0], viscos*def_grad[1], viscos*def_grad[3]}});
}

std::array<double,6>
TG2_scheme::compute_nodal_def_grad (const double& h, const double& Ux, const double& Uy, const std::array<double,2>& grad_cell_ux, const std::array<double,2>& grad_cell_uy)
{

  // def_grad = [D11, D22, D33, D12, D23, D31]

  // compute \zeta
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::sqrt( vel_x*vel_x + vel_y*vel_y );
  const double aa = h>epsilon ? 6*fluid_viscosity*abs_vel/h/yield_shear_stress : 0.;

  const double a = 3./2.;
  const double c = 65./32.;
  const double b = -(114./32.+aa);
  const double Delta = b*b-4*a*c;

  const double zeta_1 = (-b + std::sqrt(Delta))/2./a;
  const double zeta_2 = (-b - std::sqrt(Delta))/2./a;

  if ( std::abs(zeta_1 - .5)<=.5 && std::abs(zeta_2 - .5)<=.5)
  {
    std::cout << "Two valid roots, look at compute_nodal_def_grad, " << zeta_1 << " " << zeta_2 << ", STOP!" << std::endl;
    exit(1.);
  }

  const double zeta = std::abs(zeta_1 - .5)<=.5 ? zeta_1 : zeta_2;

  const auto & partial_x_ux = grad_cell_ux[0];
  const auto & partial_y_ux = grad_cell_ux[1];
  const auto   partial_z_ux = h>epsilon ? 3./(2.+zeta)*vel_x/h : 0.;

  const auto & partial_x_uy = grad_cell_uy[0];
  const auto & partial_y_uy = grad_cell_uy[1];
  const auto   partial_z_uy = h>epsilon ? 3./(2.+zeta)*vel_y/h : 0.;

  const auto partial_x_uz = 0.; // steady state simple shear flow
  const auto partial_y_uz = 0.; // steady state simple shear flow
  const auto partial_z_uz = -(grad_cell_ux[0]+grad_cell_uy[1]);
  
  return(std::array<double,6>{{partial_x_ux, partial_y_uy, partial_z_uz, .5*(partial_x_uy+partial_y_ux), .5*(partial_z_uy+partial_y_uz), .5*(partial_z_ux+partial_x_uz)}});
}


// source terms
double
TG2_scheme::h_src_formula (const double& h, const double& Ux, const double& Uy)
{ return (0.); }

double
TG2_scheme::Ux_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdx)
{
  const double bed_pressure = grav*h - surface_pressure/density;
  const double vel_x = h>epsilon*2 ? Ux/h : 0.;
  const double vel_y = h>epsilon*2 ? Uy/h : 0.;
  const double abs_vel = std::sqrt( vel_x*vel_x + vel_y*vel_y );

  const double vel_x_sign = abs_vel!=0 ? vel_x/abs_vel : 0.;

  //if (abs_vel!=0)
  //std::cout << abs_vel*vel_x << std::endl;
  //if (Ux!=0) std::cout << Ux/turbulence_coeff/(h+epsilon)/(h+epsilon)*grav*std::sqrt(Ux*Ux + Uy*Uy) << std::endl;

  const double bed_fric_contr = is_bed_friction ? (grav*abs_vel*vel_x/turbulence_coeff*0 + vel_x_sign*bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;
  return (-grav*h*dZdx - bed_fric_contr);
}

double
TG2_scheme::Uy_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdy)
{
  const double bed_pressure = grav*h - surface_pressure/density;
  const double vel_x = h>epsilon*2 ? Ux/h : 0.;
  const double vel_y = h>epsilon*2 ? Uy/h : 0.;
  const double abs_vel = std::sqrt( vel_x*vel_x + vel_y*vel_y );

  const double vel_y_sign = abs_vel!=0 ? vel_y/abs_vel : 0.;

  const double bed_fric_contr = is_bed_friction ? (grav*abs_vel*vel_y/turbulence_coeff*0 + vel_y_sign*bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;
  return (-grav*h*dZdy - bed_fric_contr);
}



