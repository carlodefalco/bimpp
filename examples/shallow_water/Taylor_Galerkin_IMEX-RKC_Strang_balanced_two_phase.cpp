#include "Taylor_Galerkin_IMEX-RKC_Strang_balanced.h"
#include <algorithm>
#include <cassert>

TG2_scheme::TG2_scheme(Q1& sol_s,
                       Q1& sol_w,
                       Q1& sold_s,
                       Q1& sold_w,
                       Q1& soldd_s, 
                       Q1& soldd_w,
                       Q1& sold_rkc_s,
                       Q1& sold_rkc_w,
                       Q1& soldd_rkc_s, 
                       Q1& soldd_rkc_w, 
                       Q1& sol_ini_rkc_s,
                       Q1& sol_ini_rkc_w,
                       Q1& incr_s,
                       Q1& incr_w,
                       Q1& incr_initial_source_s,
                       Q1& incr_initial_source_w,
                       Q1& incr_source_s,
                       Q1& incr_source_w,
                       std::vector<std::array<double,4>>& incr_anti_diff_s,
                       std::vector<std::array<double,4>>& incr_anti_diff_w,
                       Q1& stress_initial_step_s,
                       Q1& stress_initial_step_w,
                       Q1& stress_step_s,
                       Q1& stress_step_w,
                       Q1& P_plus_s,
                       Q1& P_plus_w,
                       Q1& P_minus_s,
                       Q1& P_minus_w,
                       Q1& spec_radius_nodal_s,
                       Q1& spec_radius_nodal_w,
                       Q0& sol_onehalf_s,
                       Q0& sol_onehalf_w,
                       Q1& porosity,
                       Q1& mass,
                       const ordering& oh,
                       const ordering& oUx,
                       const ordering& oUy,
                       const Q1& Z,
		                   Q1& Newton_it,
                       Q1& slope_x_node,
                       Q1& slope_y_node,
                       const Q0& slope_x,
                       const Q0& slope_y,
                       const double& DELTAT,
                       const double& h_min,
                       const bool& is_non_reflBC,
                       const bool& is_bed_friction,
                       const bool& is_stress_tensor,
                       const bool& is_erosion,
                       const double& grav,
                       const double& density,
                       const double& density_s,
                       const double& density_w,
                       const double& turbulence_coeff,
                       const double& surface_pressure, 
                       const double& bed_friction_angle_rad,
                       const double& fluid_viscosity,
                       const double& yield_shear_stress,
                       const double& erosion_coefficient)
: sol_s(sol_s), sol_w(sol_w), sold_s(sold_s), sold_w(sold_w), soldd_s(soldd_s), soldd(soldd_w), sold_rkc_s(sold_rkc_s), sold_rkc_w(sold_rkc_w), soldd_rkc_s(soldd_rkc_s), soldd_rkc_w(soldd_rkc_w), sol_ini_rkc_s(sol_ini_rkc_s), sol_ini_rkc_w(sol_ini_rkc_w), incr_s(incr_s), incr_w(incr_w), incr_initial_source_s(incr_initial_source_s), incr_initial_source_w(incr_initial_source_w), incr_source_s(incr_source_s), incr_source_w(incr_source_w), incr_anti_diff_s(incr_anti_diff_s), incr_anti_diff_w(incr_anti_diff_w), stress_initial_step_s(stress_initial_step_s), stress_initial_step_w(stress_initial_step_w), stress_step_s(stress_step_s), stress_step_w(stress_step_w), P_plus_s(P_plus_s), P_plus_w(P_plus_w), P_minus_s(P_minus_s), P_minus_w(P_minus_w), spec_radius_nodal_s(spec_radius_nodal_s), spec_radius_nodal_w(spec_radius_nodal_w), sol_onehalf_s(sol_onehalf_s), sol_onehalf_w(sol_onehalf_w), porosity(porosity), mass(mass), 
  ordh(oh), ordUx(oUx), ordUy(oUy), Z(Z), Newton_it(Newton_it), slope_x_node(slope_x_node), slope_y_node(slope_y_node), slope_x(slope_x), slope_y(slope_y), DELTAT(DELTAT), epsilon(h_min), is_non_reflBC(is_non_reflBC), is_bed_friction(is_bed_friction), is_stress_tensor(is_stress_tensor), is_erosion(is_erosion), grav(grav),
  density(density), density_s(density_s), density_w(density_w), turbulence_coeff(turbulence_coeff), surface_pressure(surface_pressure), bed_friction_angle_rad(bed_friction_angle_rad), fluid_viscosity(fluid_viscosity), yield_shear_stress(yield_shear_stress), erosion_coefficient(erosion_coefficient)
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
      hdof[ii]  = sol_s [ordh  (quadrant->gt (ii) )] + sol_w [ordh  (quadrant->gt (ii) )];
      Uxdof[ii] = sol_s [ordUx (quadrant->gt (ii) )] + sol_w [ordUx (quadrant->gt (ii) )];
      Uydof[ii] = sol_s [ordUy (quadrant->gt (ii) )] + sol_w [ordUy (quadrant->gt (ii) )];
    }
    else
    {
      hdof[ii]  = .5 * (sol_s [ordh  (quadrant->gparent (0, ii) )] +
                        sol_s [ordh  (quadrant->gparent (1, ii) )]) +
                  .5 * (sol_w [ordh  (quadrant->gparent (0, ii) )] +
                        sol_w [ordh  (quadrant->gparent (1, ii) )]);
      Uxdof[ii] = .5 * (sol_s [ordUx (quadrant->gparent (0, ii) )] +
                        sol_s [ordUx (quadrant->gparent (1, ii) )]) +
                  .5 * (sol_w [ordUx (quadrant->gparent (0, ii) )] +
                        sol_w [ordUx (quadrant->gparent (1, ii) )]);
      Uydof[ii] = .5 * (sol_s [ordUy (quadrant->gparent (0, ii) )] +
                        sol_s [ordUy (quadrant->gparent (1, ii) )]) +
                  .5 * (sol_w [ordUy (quadrant->gparent (0, ii) )] +
                        sol_w [ordUy (quadrant->gparent (1, ii) )]);
    }
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
      hdof      = sol_s   [ordh  (quadrant->gt (ii) )] + sol_w   [ordh  (quadrant->gt (ii) )];
      hdof_old  = sold_s  [ordh  (quadrant->gt (ii) )] + sold_w  [ordh  (quadrant->gt (ii) )];
      hdof_oldd = soldd_s [ordh  (quadrant->gt (ii) )] + soldd_w [ordh  (quadrant->gt (ii) )];
    }
    else
    {
      hdof      = .5 * (sol_s   [ordh  (quadrant->gparent (0, ii) )] +
                        sol_s   [ordh  (quadrant->gparent (1, ii) )]) +
                  .5 * (sol_w   [ordh  (quadrant->gparent (0, ii) )] +
                        sol_w   [ordh  (quadrant->gparent (1, ii) )]);
      hdof_old  = .5 * (sold_s  [ordh  (quadrant->gparent (0, ii) )] +
                        sold_s  [ordh  (quadrant->gparent (1, ii) )]) +
                  .5 * (sold_w  [ordh  (quadrant->gparent (0, ii) )] +
                        sold_w  [ordh  (quadrant->gparent (1, ii) )]);
      hdof_oldd = .5 * (soldd_s [ordh  (quadrant->gparent (0, ii) )] +
                        soldd_s [ordh  (quadrant->gparent (1, ii) )]) +
                  .5 * (soldd_w [ordh  (quadrant->gparent (0, ii) )] +
                        soldd_w [ordh  (quadrant->gparent (1, ii) )]);
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
  
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 
  
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii);
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;


  
  double h_cell_average_s = 0., Ux_cell_average_s = 0., Uy_cell_average_s = 0., h_cell_average_w = 0., Ux_cell_average_w = 0., Uy_cell_average_w = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    double hdof_c_s, Uxdof_c_s, Uydof_c_s, hdof_c_w, Uxdof_c_w, Uydof_c_w;

    if (! quadrant->is_hanging (ii) )
    {
      hdof_c_s    = sol_s [ordh    (quadrant->gt (ii) )];
      Uxdof_c_s   = sol_s [ordUx   (quadrant->gt (ii) )];
      Uydof_c_s   = sol_s [ordUy   (quadrant->gt (ii) )];

      hdof_c_w    = sol_w [ordh    (quadrant->gt (ii) )];
      Uxdof_c_w   = sol_w [ordUx   (quadrant->gt (ii) )];
      Uydof_c_w   = sol_w [ordUy   (quadrant->gt (ii) )];
    }
    else
    {
      hdof_c_s    = .5 * (sol_s [ordh    (quadrant->gparent (0, ii) )] +
                          sol_s [ordh    (quadrant->gparent (1, ii) )]);
      Uxdof_c_s   = .5 * (sol_s [ordUx   (quadrant->gparent (0, ii) )] +
                          sol_s [ordUx   (quadrant->gparent (1, ii) )]);
      Uydof_c_s   = .5 * (sol_s [ordUy   (quadrant->gparent (0, ii) )] +
                          sol_s [ordUy   (quadrant->gparent (1, ii) )]);

      hdof_c_w    = .5 * (sol_w [ordh    (quadrant->gparent (0, ii) )] +
                          sol_w [ordh    (quadrant->gparent (1, ii) )]);
      Uxdof_c_w   = .5 * (sol_w [ordUx   (quadrant->gparent (0, ii) )] +
                          sol_w [ordUx   (quadrant->gparent (1, ii) )]);
      Uydof_c_w   = .5 * (sol_w [ordUy   (quadrant->gparent (0, ii) )] +
                          sol_w [ordUy   (quadrant->gparent (1, ii) )]);
    }

    hdof_s [ii]   = hdof_c_s;
    Uxdof_s[ii]   = Uxdof_c_s;
    Uydof_s[ii]   = Uydof_c_s;

    hdof_w [ii]   = hdof_c_w;
    Uxdof_w[ii]   = Uxdof_c_w;
    Uydof_w[ii]   = Uydof_c_w;

    h_cell_average_s  += hdof_c_s;
    Ux_cell_average_s += Uxdof_c_s;
    Uy_cell_average_s += Uydof_c_s;

    h_cell_average_w  += hdof_c_w;
    Ux_cell_average_w += Uxdof_c_w;
    Uy_cell_average_w += Uydof_c_w;
    
    
    fluxx_h_node_s[ii]  = h_flux_formula_x   (hdof_c_s, Uxdof_c_s, Uydof_c_s);
    fluxy_h_node_s[ii]  = h_flux_formula_y   (hdof_c_s, Uxdof_c_s, Uydof_c_s);
    fluxx_Ux_node_s[ii] = Ux_flux_formula_x  (hdof_c_s, Uxdof_c_s, Uydof_c_s, hdof_c_w);
    fluxy_Ux_node_s[ii] = Ux_flux_formula_y  (hdof_c_s, Uxdof_c_s, Uydof_c_s);
    fluxx_Uy_node_s[ii] = Uy_flux_formula_x  (hdof_c_s, Uxdof_c_s, Uydof_c_s);
    fluxy_Uy_node_s[ii] = Uy_flux_formula_y  (hdof_c_s, Uxdof_c_s, Uydof_c_s, hdof_c_w);

    fluxx_h_node_w[ii]  = h_flux_formula_x   (hdof_c_w, Uxdof_c_w, Uydof_c_w);
    fluxy_h_node_w[ii]  = h_flux_formula_y   (hdof_c_w, Uxdof_c_w, Uydof_c_w);
    fluxx_Ux_node_w[ii] = Ux_flux_formula_x  (hdof_c_w, Uxdof_c_w, Uydof_c_w, hdof_c_s);
    fluxy_Ux_node_w[ii] = Ux_flux_formula_y  (hdof_c_w, Uxdof_c_w, Uydof_c_w);
    fluxx_Uy_node_w[ii] = Uy_flux_formula_x  (hdof_c_w, Uxdof_c_w, Uydof_c_w);
    fluxy_Uy_node_w[ii] = Uy_flux_formula_y  (hdof_c_w, Uxdof_c_w, Uydof_c_w, hdof_c_s);

  }
  h_cell_average_s    /= 4.;
  Ux_cell_average_s   /= 4.;
  Uy_cell_average_s   /= 4.;

  h_cell_average_w    /= 4.;
  Ux_cell_average_w   /= 4.;
  Uy_cell_average_w   /= 4.;

  
  const auto div_Fh_x_s = .5*((fluxx_h_node_s[1]-fluxx_h_node_s[0]) + (fluxx_h_node_s[3]-fluxx_h_node_s[2]));
  const auto div_Fh_y_s = .5*((fluxy_h_node_s[2]-fluxy_h_node_s[0]) + (fluxy_h_node_s[3]-fluxy_h_node_s[1]));
  const auto div_Fh_cell_s = Dy*div_Fh_x_s + Dx*div_Fh_y_s;
  
  const auto div_FUx_x_s = .5*((fluxx_Ux_node_s[1]-fluxx_Ux_node_s[0]) + (fluxx_Ux_node_s[3]-fluxx_Ux_node_s[2]));
  const auto div_FUx_y_s = .5*((fluxy_Ux_node_s[2]-fluxy_Ux_node_s[0]) + (fluxy_Ux_node_s[3]-fluxy_Ux_node_s[1]));
  const auto div_FUx_cell_s = Dy*div_FUx_x_s + Dx*div_FUx_y_s;
  
  const auto div_FUy_x_s = .5*((fluxx_Uy_node_s[1]-fluxx_Uy_node_s[0]) + (fluxx_Uy_node_s[3]-fluxx_Uy_node_s[2]));
  const auto div_FUy_y_s = .5*((fluxy_Uy_node_s[2]-fluxy_Uy_node_s[0]) + (fluxy_Uy_node_s[3]-fluxy_Uy_node_s[1]));
  const auto div_FUy_cell_s = Dy*div_FUy_x_s + Dx*div_FUy_y_s;

  //const auto h_current_s = h_cell_average_s  - (dt + dt_old)*.5*.5 *  div_Fh_cell_s /area;



  const auto div_Fh_x_w = .5*((fluxx_h_node_w[1]-fluxx_h_node_w[0]) + (fluxx_h_node_w[3]-fluxx_h_node_w[2]));
  const auto div_Fh_y_w = .5*((fluxy_h_node_w[2]-fluxy_h_node_w[0]) + (fluxy_h_node_w[3]-fluxy_h_node_w[1]));
  const auto div_Fh_cell_w = Dy*div_Fh_x_w + Dx*div_Fh_y_w;
  
  const auto div_FUx_x_w = .5*((fluxx_Ux_node_w[1]-fluxx_Ux_node_w[0]) + (fluxx_Ux_node_w[3]-fluxx_Ux_node_w[2]));
  const auto div_FUx_y_w = .5*((fluxy_Ux_node_w[2]-fluxy_Ux_node_w[0]) + (fluxy_Ux_node_w[3]-fluxy_Ux_node_w[1]));
  const auto div_FUx_cell_w = Dy*div_FUx_x_w + Dx*div_FUx_y_w;
  
  const auto div_FUy_x_w = .5*((fluxx_Uy_node_w[1]-fluxx_Uy_node_w[0]) + (fluxx_Uy_node_w[3]-fluxx_Uy_node_w[2]));
  const auto div_FUy_y_w = .5*((fluxy_Uy_node_w[2]-fluxy_Uy_node_w[0]) + (fluxy_Uy_node_w[3]-fluxy_Uy_node_w[1]));
  const auto div_FUy_cell_w = Dy*div_FUy_x_w + Dx*div_FUy_y_w;


  const double coeff_ww = -.5*density_w/density_w*grav;//*hdof_s;
  const double coeff_ws = +.5*density_w/density_w*grav;//*hdof_w;

  const double coeff_sw = +.5*density_w/density_s*grav;//*hdof_s;
  const double coeff_ss = -.5*density_w/density_s*grav;//*hdof_w;


  const double int_x1 = Dy/Dx*( (hdof_w[1] - hdof_w[0])/3. + (hdof_w[3] - hdof_w[2])/6. );
  const double int_x2 = Dy/Dx*( (hdof_w[1] - hdof_w[0])/6. + (hdof_w[3] - hdof_w[2])/3. );
  const double int_x3 = Dy/Dx*( (hdof_s[1] - hdof_s[0])/3. + (hdof_s[3] - hdof_s[2])/6. );
  const double int_x4 = Dy/Dx*( (hdof_s[1] - hdof_s[0])/6. + (hdof_s[3] - hdof_s[2])/3. );

  const double int_y1 = Dx/Dy*( (hdof_w[2] - hdof_w[0])/3. + (hdof_w[3] - hdof_w[1])/6. );
  const double int_y2 = Dx/Dy*( (hdof_w[2] - hdof_w[0])/6. + (hdof_w[3] - hdof_w[1])/3. );
  const double int_y3 = Dx/Dy*( (hdof_s[2] - hdof_s[0])/3. + (hdof_s[3] - hdof_s[1])/6. );
  const double int_y4 = Dx/Dy*( (hdof_s[2] - hdof_s[0])/6. + (hdof_s[3] - hdof_s[1])/3. );


  const double int_sw_x = .5*Dx* ( (hdof_s[0]+hdof_s[1])*int_x1 + (hdof_s[2]+hdof_s[3])*int_x2 );
  const double int_ws_x = .5*Dx* ( (hdof_w[0]+hdof_w[1])*int_x3 + (hdof_w[2]+hdof_w[3])*int_x4 );

  const double int_sw_y = .5*Dy* ( (hdof_s[0]+hdof_s[2])*int_y1 + (hdof_s[1]+hdof_s[3])*int_y2 );
  const double int_ws_y = .5*Dy* ( (hdof_w[0]+hdof_w[2])*int_y3 + (hdof_w[1]+hdof_w[3])*int_y4 );




  const auto dt_sgn = (dt + dt_old)*.5*.5;

  const auto v_h_w  = h_cell_average_w  - dt_sgn * div_Fh_cell_w/area;
  const auto v_h_s  = h_cell_average_s  - dt_sgn * div_Fh_cell_s/area;
  const auto v_Ux_w = Ux_cell_average_w - dt_sgn * ( div_FUx_cell_w + coeff_ww*int_sw_x + coeff_ws*int_ws_x )/area;
  const auto v_Ux_s = Ux_cell_average_s - dt_sgn * ( div_FUx_cell_s + coeff_sw*int_sw_x + coeff_ss*int_ws_x )/area;
  const auto v_Uy_w = Uy_cell_average_w - dt_sgn * ( div_FUy_cell_w + coeff_ww*int_sw_y + coeff_ws*int_ws_y )/area;
  const auto v_Uy_s = Uy_cell_average_s - dt_sgn * ( div_FUy_cell_s + coeff_sw*int_sw_y + coeff_ss*int_ws_y )/area;


  //const auto h_current_w = h_cell_average_w  - (dt + dt_old)*.5*.5 *  div_Fh_cell_w /area;


  /*
  sol_onehalf_s[ordh    (index_quadrant)] = h_current_s;
  sol_onehalf_s[ordUx   (index_quadrant)] = Ux_cell_average_s - (dt + dt_old)*.5*.5 * (div_FUx_cell_s/area - src_slope_formula (h_current_s, slope_x[index_quadrant]));
  sol_onehalf_s[ordUy   (index_quadrant)] = Uy_cell_average_s - (dt + dt_old)*.5*.5 * (div_FUy_cell_s/area - src_slope_formula (h_current_s, slope_y[index_quadrant]));
  
  sol_onehalf_w[ordh    (index_quadrant)] = h_current_w;
  sol_onehalf_w[ordUx   (index_quadrant)] = Ux_cell_average_w - (dt + dt_old)*.5*.5 * (div_FUx_cell_w/area - src_slope_formula (h_current_w, slope_x[index_quadrant]));
  sol_onehalf_w[ordUy   (index_quadrant)] = Uy_cell_average_w - (dt + dt_old)*.5*.5 * (div_FUy_cell_w/area - src_slope_formula (h_current_w, slope_y[index_quadrant]));
  */
  


  // solve non-linearities
  sol_onehalf_s[ordh    (index_quadrant)] = h_cell_average_s;
  sol_onehalf_s[ordUx   (index_quadrant)] = Ux_cell_average_s;
  sol_onehalf_s[ordUy   (index_quadrant)] = Uy_cell_average_s;
  
  sol_onehalf_w[ordh    (index_quadrant)] = h_cell_average_w;
  sol_onehalf_w[ordUx   (index_quadrant)] = Ux_cell_average_w;
  sol_onehalf_w[ordUy   (index_quadrant)] = Uy_cell_average_w;

  const auto & Sx_c_s     = slope_x[index_quadrant];
  const auto & Sy_c_s     = slope_y[index_quadrant];
  

  const double tolerance = 1.e-4;
  const int Nmax = 1e3;

  int count = -1;
  double error = tolerance + 1; 
  while (count++<Nmax && error>tolerance)
  {
    const auto & h_c_s  = sol_onehalf_s[ordh    (index_quadrant)];
    const auto & Ux_c_s = sol_onehalf_s[ordUx   (index_quadrant)];
    const auto & Uy_c_s = sol_onehalf_s[ordUy   (index_quadrant)];

    const auto & h_c_w  = sol_onehalf_w[ordh    (index_quadrant)];
    const auto & Ux_c_w = sol_onehalf_w[ordUx   (index_quadrant)];
    const auto & Uy_c_w = sol_onehalf_w[ordUy   (index_quadrant)];

    const auto h_c = h_c_w + h_c_s;
    const auto poro_c = h_c>epsilon ? h_c_w/h_c : 0.;

    const auto U_tot_x = Ux_c_s + Ux_c_w;
    const auto U_tot_y = Uy_c_s + Uy_c_w;

    const auto abs_mass_flux = std::sqrt(U_tot_x*U_tot_x + U_tot_y*U_tot_y);


    const auto big_A = 1. - (h_c>epsilon ? (1.-poro_c)/h_c : 0.) * dt_sgn*erosion_coefficient*abs_mass_flux + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_x)*dt_sgn*src_slope_formula (-1., Sx_c_s) + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_y)*dt_sgn*src_slope_formula (-1., Sy_c_s);
    const auto big_B = (h_c>epsilon ? poro_c/h_c : 0.)*erosion_coefficient*dt_sgn*abs_mass_flux + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_x)*dt_sgn*src_slope_formula (-1., Sx_c_s) + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_y)*dt_sgn*src_slope_formula (-1., Sy_c_s);
    const auto big_C = (h_c>epsilon ? (1.-poro_c)/h_c : 0.) * dt_sgn*erosion_coefficient*abs_mass_flux + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_x)*dt_sgn*src_slope_formula (-1., Sx_c_s) + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_y)*dt_sgn*src_slope_formula (-1., Sy_c_s);
    const auto big_D = 1. - (h_c>epsilon ? poro_c/h_c : 0.)*erosion_coefficient*dt_sgn*abs_mass_flux + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_x)*dt_sgn*src_slope_formula (-1., Sx_c_s) + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_y)*dt_sgn*src_slope_formula (-1., Sy_c_s);

    const auto F_h_w  = v_h_w  - h_c_w  + dt_sgn*h_src_formula(h_c_w, U_tot_x, U_tot_y, poro_c);
    const auto F_h_s  = v_h_s  - h_c_s  + dt_sgn*h_src_formula(h_c_w, U_tot_x, U_tot_y, 1.-poro_c);
    const auto F_Ux_w = v_Ux_w - Ux_c_w + dt_sgn*src_slope_formula (h_c_w, Sx_c_s);
    const auto F_Ux_s = v_Ux_s - Ux_c_s + dt_sgn*src_slope_formula (h_c_s, Sx_c_s);
    const auto F_Uy_w = v_Uy_w - Uy_c_w + dt_sgn*src_slope_formula (h_c_w, Sy_c_s);
    const auto F_Uy_s = v_Uy_s - Uy_c_s + dt_sgn*src_slope_formula (h_c_w, Sy_c_s);

    const auto rhs_w = F_h_w + poro_c     *erosion_coefficient*dt_sgn*signum(U_tot_x)*(F_Ux_w+F_Ux_s) + poro_c     *erosion_coefficient*dt_sgn*signum(U_tot_y)*(F_Uy_w+F_Uy_s);
    const auto rhs_s = F_h_s + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_x)*(F_Ux_w+F_Ux_s) + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_y)*(F_Uy_w+F_Uy_s);

    const auto delta_h_s = (rhs_s - C/A*rhs_w) / (big_D - big_C*big_B/big_A);
    const auto delta_h_w = (rhs_w - B*delta_h_s) / big_A;

    const auto delta_Ux_s = F_Ux_s + dt_sgn*src_slope_formula (delta_h_s, Sx_c_s);
    const auto delta_Uy_s = F_Uy_s + dt_sgn*src_slope_formula (delta_h_s, Sy_c_s);

    const auto delta_Ux_w = F_Ux_w + dt_sgn*src_slope_formula (delta_h_w, Sx_c_s);
    const auto delta_Uy_w = F_Uy_w + dt_sgn*src_slope_formula (delta_h_w, Sy_c_s);

    error = std::sqrt(delta_h_s*delta_h_s + delta_Ux_s*delta_Ux_s + delta_Uy_s*delta_Uy_s + delta_h_w*delta_h_w + delta_Ux_w*delta_Ux_w + delta_Uy_w*delta_Uy_w);

    sol_onehalf_s[ordh    (index_quadrant)] += delta_h_s ;
    sol_onehalf_s[ordUx   (index_quadrant)] += delta_Ux_s;
    sol_onehalf_s[ordUx   (index_quadrant)] += delta_Uy_s;

    sol_onehalf_w[ordh    (index_quadrant)] += delta_h_w ;
    sol_onehalf_w[ordUx   (index_quadrant)] += delta_Ux_w;
    sol_onehalf_w[ordUx   (index_quadrant)] += delta_Uy_w;

  }



}


void
TG2_scheme::solve_non_lin (const int& kk)
{

  const auto & h_old_c_s  = sol_s.       get_owned_data ()[kk  ];
  const auto & h_old_c_w  = sol_w.       get_owned_data ()[kk  ];

  const auto & Ux_old_c_s = sol_s.       get_owned_data ()[kk+1];
  const auto & Ux_old_c_w = sol_w.       get_owned_data ()[kk+1];

  const auto & Uy_old_c_s = sol_s.       get_owned_data ()[kk+2];
  const auto & Uy_old_c_w = sol_w.       get_owned_data ()[kk+2];


  const auto & Sx_c_s     = slope_x_node. get_owned_data ()[int(kk/3)];
  const auto & Sy_c_s     = slope_y_node. get_owned_data ()[int(kk/3)];

  const auto dt_sgn_ = (dt + dt_old)*.5;
  const auto dt_sgn = dt_sgn_*.5;

  const auto h_tot_old = h_old_c_w + h_old_c_s;
  const auto poro_old_c = h_tot_old>epsilon ? h_old_c_w/h_tot_old : 0.;

  const auto U_tot_old_x = Ux_old_c_s + Ux_old_c_w;
  const auto U_tot_old_y = Uy_old_c_s + Uy_old_c_w;

  const auto v_h_w  = h_old_c_w  + dt_sgn_*incr_w.get_owned_data ()[kk  ]/mass.get_owned_data ()[kk  ] + dt_sgn*h_src_formula(h_tot_old, U_tot_old_x, U_tot_old_y, poro_old_c);
  const auto v_h_s  = h_old_c_s  + dt_sgn_*incr_s.get_owned_data ()[kk  ]/mass.get_owned_data ()[kk  ] + dt_sgn*h_src_formula(h_tot_old, U_tot_old_x, U_tot_old_y, 1.-poro_old_c);
  const auto v_Ux_w = Ux_old_c_w + dt_sgn_*incr_w.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1] + dt_sgn*src_slope_formula(h_old_c_w, Sx_c_s);
  const auto v_Ux_s = Ux_old_c_s + dt_sgn_*incr_s.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1] + dt_sgn*src_slope_formula(h_old_c_s, Sx_c_s);
  const auto v_Uy_w = Uy_old_c_w + dt_sgn_*incr_w.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] + dt_sgn*src_slope_formula(h_old_c_w, Sy_c_s);
  const auto v_Uy_s = Uy_old_c_s + dt_sgn_*incr_s.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] + dt_sgn*src_slope_formula(h_old_c_s, Sy_c_s);
  

  const double tolerance = 1.e-4;
  const int Nmax = 1e3;

  int count = -1;
  double error = tolerance + 1; 
  while (count++<Nmax && error>tolerance)
  {
    const auto & h_c_s  = sol_s.get_owned_data ()[kk  ];
    const auto & Ux_c_s = sol_s.get_owned_data ()[kk+1];
    const auto & Uy_c_s = sol_s.get_owned_data ()[kk+2];

    const auto & h_c_w  = sol_w.get_owned_data ()[kk  ];
    const auto & Ux_c_w = sol_w.get_owned_data ()[kk+1];
    const auto & Uy_c_w = sol_w.get_owned_data ()[kk+2];

    const auto h_c = h_c_w + h_c_s;
    const auto poro_c = h_c>epsilon ? h_c_w/h_c : 0.;

    const auto U_tot_x = Ux_c_s + Ux_c_w;
    const auto U_tot_y = Uy_c_s + Uy_c_w;

    const auto abs_mass_flux = std::sqrt(U_tot_x*U_tot_x + U_tot_y*U_tot_y);


    const auto big_A = 1. - (h_c>epsilon ? (1.-poro_c)/h_c : 0.) * dt_sgn*erosion_coefficient*abs_mass_flux + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_x)*dt_sgn*src_slope_formula (-1., Sx_c_s) + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_y)*dt_sgn*src_slope_formula (-1., Sy_c_s);
    const auto big_B = (h_c>epsilon ? poro_c/h_c : 0.)*erosion_coefficient*dt_sgn*abs_mass_flux + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_x)*dt_sgn*src_slope_formula (-1., Sx_c_s) + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_y)*dt_sgn*src_slope_formula (-1., Sy_c_s);
    const auto big_C = (h_c>epsilon ? (1.-poro_c)/h_c : 0.) * dt_sgn*erosion_coefficient*abs_mass_flux + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_x)*dt_sgn*src_slope_formula (-1., Sx_c_s) + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_y)*dt_sgn*src_slope_formula (-1., Sy_c_s);
    const auto big_D = 1. - (h_c>epsilon ? poro_c/h_c : 0.)*erosion_coefficient*dt_sgn*abs_mass_flux + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_x)*dt_sgn*src_slope_formula (-1., Sx_c_s) + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_y)*dt_sgn*src_slope_formula (-1., Sy_c_s);

    const auto F_h_w  = v_h_w  - h_c_w  + dt_sgn*h_src_formula(h_c_w, U_tot_x, U_tot_y, poro_c);
    const auto F_h_s  = v_h_s  - h_c_s  + dt_sgn*h_src_formula(h_c_w, U_tot_x, U_tot_y, 1.-poro_c);
    const auto F_Ux_w = v_Ux_w - Ux_c_w + dt_sgn*src_slope_formula (h_c_w, Sx_c_s);
    const auto F_Ux_s = v_Ux_s - Ux_c_s + dt_sgn*src_slope_formula (h_c_s, Sx_c_s);
    const auto F_Uy_w = v_Uy_w - Uy_c_w + dt_sgn*src_slope_formula (h_c_w, Sy_c_s);
    const auto F_Uy_s = v_Uy_s - Uy_c_s + dt_sgn*src_slope_formula (h_c_w, Sy_c_s);

    const auto rhs_w = F_h_w + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_x)*(F_Ux_w+F_Ux_s) + poro_c*erosion_coefficient*dt_sgn*signum(U_tot_y)*(F_Uy_w+F_Uy_s);
    const auto rhs_s = F_h_s + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_x)*(F_Ux_w+F_Ux_s) + (1.-poro_c)*erosion_coefficient*dt_sgn*signum(U_tot_y)*(F_Uy_w+F_Uy_s);

    const auto delta_h_s = (rhs_s - C/A*rhs_w) / (big_D - big_C*big_B/big_A);
    const auto delta_h_w = (rhs_w - B*delta_h_s) / big_A;

    const auto delta_Ux_s = F_Ux_s + dt_sgn*src_slope_formula (delta_h_s, Sx_c_s);
    const auto delta_Uy_s = F_Uy_s + dt_sgn*src_slope_formula (delta_h_s, Sy_c_s);

    const auto delta_Ux_w = F_Ux_w + dt_sgn*src_slope_formula (delta_h_w, Sx_c_s);
    const auto delta_Uy_w = F_Uy_w + dt_sgn*src_slope_formula (delta_h_w, Sy_c_s);

    error = std::sqrt(delta_h_s*delta_h_s + delta_Ux_s*delta_Ux_s + delta_Uy_s*delta_Uy_s + delta_h_w*delta_h_w + delta_Ux_w*delta_Ux_w + delta_Uy_w*delta_Uy_w);


    sol_s.get_owned_data ()[kk  ] += delta_h_s; 
    sol_w.get_owned_data ()[kk  ] += delta_h_w;

    sol_s.get_owned_data ()[kk+1] += delta_Ux_s; 
    sol_w.get_owned_data ()[kk+1] += delta_Ux_w;

    sol_s.get_owned_data ()[kk+2] += delta_Uy_s; 
    sol_w.get_owned_data ()[kk+2] += delta_Uy_w;

  }


}


void
TG2_scheme::compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant)
{

  // look at tmesh.h
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 

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
      hdof_s [ii]   = sol_s [ordh    (quadrant->gt (ii))];
      Uxdof_s[ii]   = sol_s [ordUx   (quadrant->gt (ii))];
      Uydof_s[ii]   = sol_s [ordUy   (quadrant->gt (ii))];

      hdof_w [ii]   = sol_w [ordh    (quadrant->gt (ii))];
      Uxdof_w[ii]   = sol_w [ordUx   (quadrant->gt (ii))];
      Uydof_w[ii]   = sol_w [ordUy   (quadrant->gt (ii))];
      
      isdof_or_hanging[ii] = 1.;
    } else {
      hdof_s [ii]   = .5 * (sol_s [ordh  (quadrant->gparent(0,ii))] +
                            sol_s [ordh  (quadrant->gparent(1,ii))]);
      Uxdof_s[ii]   = .5 * (sol_s [ordUx (quadrant->gparent(0,ii))] +
                            sol_s [ordUx (quadrant->gparent(1,ii))]);
      Uydof_s[ii]   = .5 * (sol_s [ordUy (quadrant->gparent(0,ii))] +
                            sol_s [ordUy (quadrant->gparent(1,ii))]);

      hdof_w [ii]   = .5 * (sol_w [ordh  (quadrant->gparent(0,ii))] +
                            sol_w [ordh  (quadrant->gparent(1,ii))]);
      Uxdof_w[ii]   = .5 * (sol_w [ordUx (quadrant->gparent(0,ii))] +
                            sol_w [ordUx (quadrant->gparent(1,ii))]);
      Uydof_w[ii]   = .5 * (sol_w [ordUy (quadrant->gparent(0,ii))] +
                            sol_w [ordUy (quadrant->gparent(1,ii))]);
      
      isdof_or_hanging[ii] = .5;
    }
    hdof [ii]   = hdof_s [ii] + hdof_w [ii];
    Uxdof[ii]   = Uxdof_s[ii] + Uxdof_w[ii];
    Uydof[ii]   = Uydof_s[ii] + Uydof_w[ii];

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

  //const std::array<double,4> eta_vec = {hdof[0]+Z_node[0], hdof[1]+Z_node[1], hdof[2]+Z_node[2], hdof[3]+Z_node[3]};
  //grad_cell_eta  = {.5 * ( (eta_vec[3] - eta_vec[2]) + (eta_vec[1] - eta_vec[0]) ), .5 * ( (eta_vec[2] - eta_vec[0]) + (eta_vec[3] - eta_vec[1]) )};
  
  grad_cell_eta_s  = {.5 * ( (hdof_s   [3] - hdof_s   [2]) + (hdof_s   [1] - hdof_s   [0]) ), .5 * ( (hdof_s   [2] - hdof_s   [0]) + (hdof_s   [3] - hdof_s   [1]) )};
  grad_cell_Ux_s   = {.5 * ( (Uxdof_s  [3] - Uxdof_s  [2]) + (Uxdof_s  [1] - Uxdof_s  [0]) ), .5 * ( (Uxdof_s  [2] - Uxdof_s  [0]) + (Uxdof_s  [3] - Uxdof_s  [1]) )};
  grad_cell_Uy_s   = {.5 * ( (Uydof_s  [3] - Uydof_s  [2]) + (Uydof_s  [1] - Uydof_s  [0]) ), .5 * ( (Uydof_s  [2] - Uydof_s  [0]) + (Uydof_s  [3] - Uydof_s  [1]) )};

  grad_cell_eta_w  = {.5 * ( (hdof_w   [3] - hdof_w   [2]) + (hdof_w   [1] - hdof_w   [0]) ), .5 * ( (hdof_w   [2] - hdof_w   [0]) + (hdof_w   [3] - hdof_w   [1]) )};
  grad_cell_Ux_w   = {.5 * ( (Uxdof_w  [3] - Uxdof_w  [2]) + (Uxdof_w  [1] - Uxdof_w  [0]) ), .5 * ( (Uxdof_w  [2] - Uxdof_w  [0]) + (Uxdof_w  [3] - Uxdof_w  [1]) )};
  grad_cell_Uy_w   = {.5 * ( (Uydof_w  [3] - Uydof_w  [2]) + (Uydof_w  [1] - Uydof_w  [0]) ), .5 * ( (Uydof_w  [2] - Uydof_w  [0]) + (Uydof_w  [3] - Uydof_w  [1]) )};



  //std::cout << grad_cell_eta[0] << " " << grad_cell_eta[1] << std::endl;

  const double & h_cell_s    = sol_onehalf_s[ordh    (index_quadrant)];
  const double & Ux_cell_s   = sol_onehalf_s[ordUx   (index_quadrant)];
  const double & Uy_cell_s   = sol_onehalf_s[ordUy   (index_quadrant)];

  const double & h_cell_w    = sol_onehalf_w[ordh    (index_quadrant)];
  const double & Ux_cell_w   = sol_onehalf_w[ordUx   (index_quadrant)];
  const double & Uy_cell_w   = sol_onehalf_w[ordUy   (index_quadrant)];

  std::array<double,4> contr_x_1 = {0., 0., 0., 0.}, contr_x_2 = {0., 0., 0., 0.}, contr_y_1 = {0., 0., 0., 0.}, contr_y_2 = {0., 0., 0., 0.};



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
      }

      const auto & index_quadrant_nei = quadrant_nei->get_forest_quad_idx (); 

      for (int jEdge = 0; jEdge < 4; ++jEdge) { // cycle neigh edges 

        const auto j_1 = bimpp_to_rev_ord[jEdge];
        const auto j_2 = bimpp_to_rev_ord[(jEdge+1)%4];

        const auto edge_length_nei = std::sqrt(std::pow((Xn[j_1]-Xn[j_2]),2.) + std::pow((Yn[j_1]-Yn[j_2]),2.));
        const std::array<double,2> outward_normal_edge_nei = {(-Yn[j_1]+Yn[j_2])/edge_length_nei, ( Xn[j_1]-Xn[j_2])/edge_length_nei};  
        const bool check_orthogonality = std::inner_product(outward_normal_edge_nei.begin(), outward_normal_edge_nei.end(), outward_normal_edge.begin(), 0.) == -1;


        if ( (((xn[i_1] == Xn[j_1] && yn[i_1] == Yn[j_1]) || 
               (xn[i_2] == Xn[j_1] && yn[i_2] == Yn[j_1]))||
              ((xn[i_1] == Xn[j_2] && yn[i_1] == Yn[j_2]) ||
               (xn[i_2] == Xn[j_2] && yn[i_2] == Yn[j_2]))) && check_orthogonality && index_quadrant!=index_quadrant_nei )
        {
          is_boundary_edge = false;

          // scrivere qui la somma dei contributi per i termini non-cons.!
          auto h_cell_w_nei = sol_onehalf_w[ordh(index_quadrant_nei)];
          auto h_cell_s_nei = sol_onehalf_s[ordh(index_quadrant_nei)];


          contr_x_1[i_1] += .5*signum(outward_normal_edge[0])*h_cell_s*(h_cell_w_nei - h_cell_w)*isdof_or_hanging[i_1];
          contr_x_2[i_1] += .5*signum(outward_normal_edge[0])*h_cell_w*(h_cell_s_nei - h_cell_s)*isdof_or_hanging[i_1];

          contr_x_1[i_2] += .5*signum(outward_normal_edge[0])*h_cell_s*(h_cell_w_nei - h_cell_w)*isdof_or_hanging[i_2];
          contr_x_2[i_2] += .5*signum(outward_normal_edge[0])*h_cell_w*(h_cell_s_nei - h_cell_s)*isdof_or_hanging[i_2];


          contr_y_1[i_1] += .5*signum(outward_normal_edge[1])*h_cell_s*(h_cell_w_nei - h_cell_w)*isdof_or_hanging[i_1];
          contr_y_2[i_1] += .5*signum(outward_normal_edge[1])*h_cell_w*(h_cell_s_nei - h_cell_s)*isdof_or_hanging[i_1];

          contr_y_1[i_2] += .5*signum(outward_normal_edge[1])*h_cell_s*(h_cell_w_nei - h_cell_w)*isdof_or_hanging[i_2];
          contr_y_2[i_2] += .5*signum(outward_normal_edge[1])*h_cell_w*(h_cell_s_nei - h_cell_s)*isdof_or_hanging[i_2];



          //break; // this just goes outside the jEdge cycle 
        }
      }

    }

    if (is_boundary_edge) // set boundary conditions
    { 
      
      auto h_cell_nei_s  = h_cell_s;
      auto Ux_cell_nei_s = Ux_cell_s;
      auto Uy_cell_nei_s = Uy_cell_s;

      auto h_cell_nei_w  = h_cell_w;
      auto Ux_cell_nei_w = Ux_cell_w;
      auto Uy_cell_nei_w = Uy_cell_w;



      Ux_cell_nei_s -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell_s + outward_normal_edge[1]*Uy_cell_s)*outward_normal_edge[0];
      Uy_cell_nei_s -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell_s + outward_normal_edge[1]*Uy_cell_s)*outward_normal_edge[1];

      const auto speed_s     = h_cell_s    >epsilon ? std::abs((Ux_cell_s    /h_cell_s    )*outward_normal_edge[0] + (Uy_cell_s    /h_cell_s    )*outward_normal_edge[1]) + std::sqrt(grav*h_cell_s    ) : 0.;
      const auto speed_nei_s = h_cell_nei_s>epsilon ? std::abs((Ux_cell_nei_s/h_cell_nei_s)*outward_normal_edge[0] + (Uy_cell_nei_s/h_cell_nei_s)*outward_normal_edge[1]) + std::sqrt(grav*h_cell_nei_s) : 0.;

      const auto smax_s = std::max(speed_s, speed_nei_s); 

      const auto flux_int_h_s  = .5*((h_flux_formula_x (h_cell_s, Ux_cell_s, Uy_cell_s)          +h_flux_formula_x (h_cell_nei_s, Ux_cell_nei_s, Uy_cell_nei_s))*outward_normal_edge[0]               + (h_flux_formula_y (h_cell_s, Ux_cell_s, Uy_cell_s)          +h_flux_formula_y (h_cell_nei_s, Ux_cell_nei_s, Uy_cell_nei_s))              *outward_normal_edge[1]) - .5*smax_s*(h_cell_nei_s -h_cell_s );
      const auto flux_int_Ux_s = .5*((Ux_flux_formula_x(h_cell_s, Ux_cell_s, Uy_cell_s, h_cell_w)+Ux_flux_formula_x(h_cell_nei_s, Ux_cell_nei_s, Uy_cell_nei_s, h_cell_nei_w))*outward_normal_edge[0] + (Ux_flux_formula_y(h_cell_s, Ux_cell_s, Uy_cell_s)          +Ux_flux_formula_y(h_cell_nei_s, Ux_cell_nei_s, Uy_cell_nei_s))              *outward_normal_edge[1]) - .5*smax_s*(Ux_cell_nei_s-Ux_cell_s);
      const auto flux_int_Uy_s = .5*((Uy_flux_formula_x(h_cell_s, Ux_cell_s, Uy_cell_s)          +Uy_flux_formula_x(h_cell_nei_s, Ux_cell_nei_s, Uy_cell_nei_s))*outward_normal_edge[0]               + (Uy_flux_formula_y(h_cell_s, Ux_cell_s, Uy_cell_s, h_cell_w)+Uy_flux_formula_y(h_cell_nei_s, Ux_cell_nei_s, Uy_cell_nei_s, h_cell_nei_w))*outward_normal_edge[1]) - .5*smax_s*(Uy_cell_nei_s-Uy_cell_s);




      Ux_cell_nei_w -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell_w + outward_normal_edge[1]*Uy_cell_w)*outward_normal_edge[0];
      Uy_cell_nei_w -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell_w + outward_normal_edge[1]*Uy_cell_w)*outward_normal_edge[1];

      const auto speed_w     = h_cell_w    >epsilon ? std::abs((Ux_cell_w    /h_cell_w    )*outward_normal_edge[0] + (Uy_cell_w    /h_cell_w    )*outward_normal_edge[1]) + std::sqrt(grav*h_cell_w    ) : 0.;
      const auto speed_nei_w = h_cell_nei_w>epsilon ? std::abs((Ux_cell_nei_w/h_cell_nei_w)*outward_normal_edge[0] + (Uy_cell_nei_w/h_cell_nei_w)*outward_normal_edge[1]) + std::sqrt(grav*h_cell_nei_w) : 0.;

      const auto smax_w = std::max(speed_w, speed_nei_w); 

      const auto flux_int_h_w  = .5*((h_flux_formula_x (h_cell_w, Ux_cell_w, Uy_cell_w)          +h_flux_formula_x (h_cell_nei_w, Ux_cell_nei_w, Uy_cell_nei_w))*outward_normal_edge[0]               + (h_flux_formula_y (h_cell_w, Ux_cell_w, Uy_cell_w)          +h_flux_formula_y (h_cell_nei_w, Ux_cell_nei_w, Uy_cell_nei_w))              *outward_normal_edge[1]) - .5*smax_w*(h_cell_nei_w -h_cell_w );
      const auto flux_int_Ux_w = .5*((Ux_flux_formula_x(h_cell_w, Ux_cell_w, Uy_cell_w, h_cell_s)+Ux_flux_formula_x(h_cell_nei_w, Ux_cell_nei_w, Uy_cell_nei_w, h_cell_nei_s))*outward_normal_edge[0] + (Ux_flux_formula_y(h_cell_w, Ux_cell_w, Uy_cell_w)          +Ux_flux_formula_y(h_cell_nei_w, Ux_cell_nei_w, Uy_cell_nei_w))              *outward_normal_edge[1]) - .5*smax_w*(Ux_cell_nei_w-Ux_cell_w);
      const auto flux_int_Uy_w = .5*((Uy_flux_formula_x(h_cell_w, Ux_cell_w, Uy_cell_w)          +Uy_flux_formula_x(h_cell_nei_w, Ux_cell_nei_w, Uy_cell_nei_w))*outward_normal_edge[0]               + (Uy_flux_formula_y(h_cell_w, Ux_cell_w, Uy_cell_w, h_cell_s)+Uy_flux_formula_y(h_cell_nei_w, Ux_cell_nei_w, Uy_cell_nei_w, h_cell_nei_s))*outward_normal_edge[1]) - .5*smax_w*(Uy_cell_nei_w-Uy_cell_w);


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
        incr_s[ordh   (quadrant->gt (i_1))] += -edge_length*flux_int_h_s *.5;
        incr_s[ordUx  (quadrant->gt (i_1))] += -edge_length*flux_int_Ux_s*.5;
        incr_s[ordUy  (quadrant->gt (i_1))] += -edge_length*flux_int_Uy_s*.5;
      }
      else
      {
        incr_s [ordh  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_h_s *.5*.5;
        incr_s [ordh  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_h_s *.5*.5;
      
        incr_s [ordUx (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Ux_s*.5*.5;
        incr_s [ordUx (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Ux_s*.5*.5;
      
        incr_s [ordUy (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uy_s*.5*.5;
        incr_s [ordUy (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uy_s*.5*.5;
      }

      if (! quadrant->is_hanging (i_2))
      {
        incr_s[ordh   (quadrant->gt (i_2))] += -edge_length*flux_int_h_s *.5;
        incr_s[ordUx  (quadrant->gt (i_2))] += -edge_length*flux_int_Ux_s*.5;
        incr_s[ordUy  (quadrant->gt (i_2))] += -edge_length*flux_int_Uy_s*.5;
      }
      else
      {
        // il secondo .5 è per hanging nodes
        incr_s [ordh  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_h_s *.5*.5;
        incr_s [ordh  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_h_s *.5*.5;
      
        incr_s [ordUx (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Ux_s*.5*.5;
        incr_s [ordUx (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Ux_s*.5*.5;
      
        incr_s [ordUy (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uy_s*.5*.5;
        incr_s [ordUy (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uy_s*.5*.5;
      }



      // .5 is the base function evaluated in the middle, mid-point intergration
      if (! quadrant->is_hanging (i_1))
      {
        incr_w[ordh   (quadrant->gt (i_1))] += -edge_length*flux_int_h_w *.5;
        incr_w[ordUx  (quadrant->gt (i_1))] += -edge_length*flux_int_Ux_w*.5;
        incr_w[ordUy  (quadrant->gt (i_1))] += -edge_length*flux_int_Uy_w*.5;
      }
      else
      {
        incr_w [ordh  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_h_w *.5*.5;
        incr_w [ordh  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_h_w *.5*.5;
      
        incr_w [ordUx (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Ux_w*.5*.5;
        incr_w [ordUx (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Ux_w*.5*.5;
      
        incr_w [ordUy (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uy_w*.5*.5;
        incr_w [ordUy (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uy_w*.5*.5;
      }

      if (! quadrant->is_hanging (i_2))
      {
        incr_w[ordh   (quadrant->gt (i_2))] += -edge_length*flux_int_h_w *.5;
        incr_w[ordUx  (quadrant->gt (i_2))] += -edge_length*flux_int_Ux_w*.5;
        incr_w[ordUy  (quadrant->gt (i_2))] += -edge_length*flux_int_Uy_w*.5;
      }
      else
      {
        // il secondo .5 è per hanging nodes
        incr_w [ordh  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_h_w *.5*.5;
        incr_w [ordh  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_h_w *.5*.5;
      
        incr_w [ordUx (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Ux_w*.5*.5;
        incr_w [ordUx (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Ux_w*.5*.5;
      
        incr_w [ordUy (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uy_w*.5*.5;
        incr_w [ordUy (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uy_w*.5*.5;
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
  

  const auto diff_term_h_x_s  = grad_cell_eta_s[0] * vel_rusanov_cell_y*.5;
  const auto diff_term_h_y_s  = grad_cell_eta_s[1] * vel_rusanov_cell_x*.5;

  const auto diff_term_Ux_x_s = grad_cell_Ux_s [0] * vel_rusanov_cell_y*.5;
  const auto diff_term_Ux_y_s = grad_cell_Ux_s [1] * vel_rusanov_cell_x*.5;

  const auto diff_term_Uy_x_s = grad_cell_Uy_s [0] * vel_rusanov_cell_y*.5;
  const auto diff_term_Uy_y_s = grad_cell_Uy_s [1] * vel_rusanov_cell_x*.5;


  const auto F_star_h_x_s  = h_flux_formula_x(h_cell_s, Ux_cell_s, Uy_cell_s)            - diff_term_h_x_s;
  const auto F_star_h_y_s  = h_flux_formula_y(h_cell_s, Ux_cell_s, Uy_cell_s)            - diff_term_h_y_s;

  const auto F_star_Ux_x_s = Ux_flux_formula_x(h_cell_s, Ux_cell_s, Uy_cell_s, h_cell_w) - diff_term_Ux_x_s;
  const auto F_star_Ux_y_s = Ux_flux_formula_y(h_cell_s, Ux_cell_s, Uy_cell_s)           - diff_term_Ux_y_s;

  const auto F_star_Uy_x_s = Uy_flux_formula_x(h_cell_s, Ux_cell_s, Uy_cell_s)           - diff_term_Uy_x_s;
  const auto F_star_Uy_y_s = Uy_flux_formula_y(h_cell_s, Ux_cell_s, Uy_cell_s, h_cell_w) - diff_term_Uy_y_s;




  const auto diff_term_h_x_w  = grad_cell_eta_w[0] * vel_rusanov_cell_y*.5;
  const auto diff_term_h_y_w  = grad_cell_eta_w[1] * vel_rusanov_cell_x*.5;

  const auto diff_term_Ux_x_w = grad_cell_Ux_w [0] * vel_rusanov_cell_y*.5;
  const auto diff_term_Ux_y_w = grad_cell_Ux_w [1] * vel_rusanov_cell_x*.5;

  const auto diff_term_Uy_x_w = grad_cell_Uy_w [0] * vel_rusanov_cell_y*.5;
  const auto diff_term_Uy_y_w = grad_cell_Uy_w [1] * vel_rusanov_cell_x*.5;


  const auto F_star_h_x_w  = h_flux_formula_x(h_cell_w, Ux_cell_w, Uy_cell_w)            - diff_term_h_x_w;
  const auto F_star_h_y_w  = h_flux_formula_y(h_cell_w, Ux_cell_w, Uy_cell_w)            - diff_term_h_y_w;

  const auto F_star_Ux_x_w = Ux_flux_formula_x(h_cell_w, Ux_cell_w, Uy_cell_w, h_cell_s) - diff_term_Ux_x_w;
  const auto F_star_Ux_y_w = Ux_flux_formula_y(h_cell_w, Ux_cell_w, Uy_cell_w)           - diff_term_Ux_y_w;

  const auto F_star_Uy_x_w = Uy_flux_formula_x(h_cell_w, Ux_cell_w, Uy_cell_w)           - diff_term_Uy_x_w;
  const auto F_star_Uy_y_w = Uy_flux_formula_y(h_cell_w, Ux_cell_w, Uy_cell_w, h_cell_s) - diff_term_Uy_y_w;



  const double coeff_ww = -.5*density_w/density_w*grav;//*hdof_s;
  const double coeff_ws = +.5*density_w/density_w*grav;//*hdof_w;

  const double coeff_sw = +.5*density_w/density_s*grav;//*hdof_s;
  const double coeff_ss = -.5*density_w/density_s*grav;//*hdof_w;



  for (int ii = 0; ii < 4; ++ii){

    const auto h_s_    = der_coeffs_x[ii]*F_star_h_x_s +der_coeffs_y[ii]*F_star_h_y_s;
    const auto Ux_s_   = der_coeffs_x[ii]*F_star_Ux_x_s+der_coeffs_y[ii]*F_star_Ux_y_s + coeff_sw*Dy*contr_x_1[ii] + coeff_ss*Dy*contr_x_2[ii];
    const auto Uy_s_   = der_coeffs_x[ii]*F_star_Uy_x_s+der_coeffs_y[ii]*F_star_Uy_y_s + coeff_sw*Dx*contr_y_1[ii] + coeff_ss*Dx*contr_y_2[ii];

    const auto h_al_s  = der_coeffs_x[ii]*diff_term_h_x_s  + der_coeffs_y[ii]*diff_term_h_y_s; 
    const auto Ux_al_s = der_coeffs_x[ii]*diff_term_Ux_x_s + der_coeffs_y[ii]*diff_term_Ux_y_s;
    const auto Uy_al_s = der_coeffs_x[ii]*diff_term_Uy_x_s + der_coeffs_y[ii]*diff_term_Uy_y_s;


    const auto h_w_    = der_coeffs_x[ii]*F_star_h_x_w +der_coeffs_y[ii]*F_star_h_y_w;
    const auto Ux_w_   = der_coeffs_x[ii]*F_star_Ux_x_w+der_coeffs_y[ii]*F_star_Ux_y_w + coeff_ww*Dy*contr_x_1[ii] + coeff_ws*Dy*contr_x_2[ii];
    const auto Uy_w_   = der_coeffs_x[ii]*F_star_Uy_x_w+der_coeffs_y[ii]*F_star_Uy_y_w + coeff_ww*Dx*contr_y_1[ii] + coeff_ws*Dx*contr_y_2[ii];

    const auto h_al_w  = der_coeffs_x[ii]*diff_term_h_x_w  + der_coeffs_y[ii]*diff_term_h_y_w; 
    const auto Ux_al_w = der_coeffs_x[ii]*diff_term_Ux_x_w + der_coeffs_y[ii]*diff_term_Ux_y_w;
    const auto Uy_al_w = der_coeffs_x[ii]*diff_term_Uy_x_w + der_coeffs_y[ii]*diff_term_Uy_y_w;


    incr_anti_diff_s[ordh (index_quadrant)][ii] = h_al_s;
    incr_anti_diff_s[ordUx(index_quadrant)][ii] = Ux_al_s;
    incr_anti_diff_s[ordUy(index_quadrant)][ii] = Uy_al_s;

    incr_anti_diff_w[ordh (index_quadrant)][ii] = h_al_w;
    incr_anti_diff_w[ordUx(index_quadrant)][ii] = Ux_al_w;
    incr_anti_diff_w[ordUy(index_quadrant)][ii] = Uy_al_w;


    if (! quadrant->is_hanging (ii)){

      incr_s [ordh  (quadrant->gt (ii))] += h_s_;
      incr_s [ordUx (quadrant->gt (ii))] += Ux_s_;
      incr_s [ordUy (quadrant->gt (ii))] += Uy_s_;

      P_plus_s [ordh  (quadrant->gt (ii))] += std::max(0., h_al_s );
      P_plus_s [ordUx (quadrant->gt (ii))] += std::max(0., Ux_al_s);
      P_plus_s [ordUy (quadrant->gt (ii))] += std::max(0., Uy_al_s);

      P_minus_s [ordh  (quadrant->gt (ii))] += std::min(0., h_al_s );
      P_minus_s [ordUx (quadrant->gt (ii))] += std::min(0., Ux_al_s);
      P_minus_s [ordUy (quadrant->gt (ii))] += std::min(0., Uy_al_s);




      incr_w [ordh  (quadrant->gt (ii))] += h_w_;
      incr_w [ordUx (quadrant->gt (ii))] += Ux_w_;
      incr_w [ordUy (quadrant->gt (ii))] += Uy_w_;

      P_plus_w [ordh  (quadrant->gt (ii))] += std::max(0., h_al_w );
      P_plus_w [ordUx (quadrant->gt (ii))] += std::max(0., Ux_al_w);
      P_plus_w [ordUy (quadrant->gt (ii))] += std::max(0., Uy_al_w);

      P_minus_w [ordh  (quadrant->gt (ii))] += std::min(0., h_al_w );
      P_minus_w [ordUx (quadrant->gt (ii))] += std::min(0., Ux_al_w);
      P_minus_w [ordUy (quadrant->gt (ii))] += std::min(0., Uy_al_w);
      
      
    } else {


      incr_s [ordh  (quadrant->gparent(0,ii))] += h_s_;
      incr_s [ordh  (quadrant->gparent(1,ii))] += h_s_;
      
      incr_s [ordUx (quadrant->gparent(0,ii))] += Ux_s_;
      incr_s [ordUx (quadrant->gparent(1,ii))] += Ux_s_;
      
      incr_s [ordUy (quadrant->gparent(0,ii))] += Uy_s_;
      incr_s [ordUy (quadrant->gparent(1,ii))] += Uy_s_;



      P_plus_s [ordh  (quadrant->gparent(0,ii))] += std::max(0., h_al_s);
      P_plus_s [ordh  (quadrant->gparent(1,ii))] += std::max(0., h_al_s);
      
      P_plus_s [ordUx (quadrant->gparent(0,ii))] += std::max(0., Ux_al_s);
      P_plus_s [ordUx (quadrant->gparent(1,ii))] += std::max(0., Ux_al_s);
      
      P_plus_s [ordUy (quadrant->gparent(0,ii))] += std::max(0., Uy_al_s);
      P_plus_s [ordUy (quadrant->gparent(1,ii))] += std::max(0., Uy_al_s);



      P_minus_s [ordh  (quadrant->gparent(0,ii))] += std::min(0., h_al_s);
      P_minus_s [ordh  (quadrant->gparent(1,ii))] += std::min(0., h_al_s);
      
      P_minus_s [ordUx (quadrant->gparent(0,ii))] += std::min(0., Ux_al_s);
      P_minus_s [ordUx (quadrant->gparent(1,ii))] += std::min(0., Ux_al_s);
      
      P_minus_s [ordUy (quadrant->gparent(0,ii))] += std::min(0., Uy_al_s);
      P_minus_s [ordUy (quadrant->gparent(1,ii))] += std::min(0., Uy_al_s);







      incr_w [ordh  (quadrant->gparent(0,ii))] += h_w_;
      incr_w [ordh  (quadrant->gparent(1,ii))] += h_w_;
      
      incr_w [ordUx (quadrant->gparent(0,ii))] += Ux_w_;
      incr_w [ordUx (quadrant->gparent(1,ii))] += Ux_w_;
      
      incr_w [ordUy (quadrant->gparent(0,ii))] += Uy_w_;
      incr_w [ordUy (quadrant->gparent(1,ii))] += Uy_w_;



      P_plus_w [ordh  (quadrant->gparent(0,ii))] += std::max(0., h_al_w);
      P_plus_w [ordh  (quadrant->gparent(1,ii))] += std::max(0., h_al_w);
      
      P_plus_w [ordUx (quadrant->gparent(0,ii))] += std::max(0., Ux_al_w);
      P_plus_w [ordUx (quadrant->gparent(1,ii))] += std::max(0., Ux_al_w);
      
      P_plus_w [ordUy (quadrant->gparent(0,ii))] += std::max(0., Uy_al_w);
      P_plus_w [ordUy (quadrant->gparent(1,ii))] += std::max(0., Uy_al_w);



      P_minus_w [ordh  (quadrant->gparent(0,ii))] += std::min(0., h_al_w);
      P_minus_w [ordh  (quadrant->gparent(1,ii))] += std::min(0., h_al_w);
      
      P_minus_w [ordUx (quadrant->gparent(0,ii))] += std::min(0., Ux_al_w);
      P_minus_w [ordUx (quadrant->gparent(1,ii))] += std::min(0., Ux_al_w);
      
      P_minus_w [ordUy (quadrant->gparent(0,ii))] += std::min(0., Uy_al_w);
      P_minus_w [ordUy (quadrant->gparent(1,ii))] += std::min(0., Uy_al_w);
      
    }


  }



}


void
TG2_scheme::loop_step (const int& kk, const bool& isInitial) //
{

  const auto & h_c_s  = soldd_rkc_s.get_owned_data ()[kk  ];
  const auto & Ux_c_s = soldd_rkc_s.get_owned_data ()[kk+1];
  const auto & Uy_c_s = soldd_rkc_s.get_owned_data ()[kk+2]; 

  const auto & h_c_w  = soldd_rkc_w.get_owned_data ()[kk  ];
  const auto & Ux_c_w = soldd_rkc_w.get_owned_data ()[kk+1];
  const auto & Uy_c_w = soldd_rkc_w.get_owned_data ()[kk+2]; 

  const auto & poro_c = porosity.get_owned_data ()[int(kk/3)];


  const auto h_c = h_c_s + h_c_w;
    

  const auto h_s_s_  = 0.;
  const auto Ux_s_s_ = Ux_src_formula_s(h_c, Ux_c_w, Ux_c_s, Uy_c_w, Uy_c_s, poro_c); 
  const auto Uy_s_s_ = Uy_src_formula_s(h_c, Ux_c_w, Ux_c_s, Uy_c_w, Uy_c_s, poro_c); 

  const auto h_s_w_  = 0.;
  const auto Ux_s_w_ = Ux_src_formula_w(h_c, Ux_c_w, Ux_c_s, Uy_c_w, Uy_c_s, poro_c); 
  const auto Uy_s_w_ = Uy_src_formula_w(h_c, Ux_c_w, Ux_c_s, Uy_c_w, Uy_c_s, poro_c); 

  
  if (isInitial)
  {
    incr_initial_source_s.get_owned_data ()[kk  ] = h_s_s_; 
    incr_initial_source_s.get_owned_data ()[kk+1] = Ux_s_s_;
    incr_initial_source_s.get_owned_data ()[kk+2] = Uy_s_s_;   

    incr_initial_source_w.get_owned_data ()[kk  ] = h_s_w_; 
    incr_initial_source_w.get_owned_data ()[kk+1] = Ux_s_w_;
    incr_initial_source_w.get_owned_data ()[kk+2] = Uy_s_w_;    
  }
  else
  {
    incr_source_s.get_owned_data ()[kk  ] = h_s_s_; 
    incr_source_s.get_owned_data ()[kk+1] = Ux_s_s_;
    incr_source_s.get_owned_data ()[kk+2] = Uy_s_s_;

    incr_source_w.get_owned_data ()[kk  ] = h_s_w_; 
    incr_source_w.get_owned_data ()[kk+1] = Ux_s_w_;
    incr_source_w.get_owned_data ()[kk+2] = Uy_s_w_;
  }

}

/*
void
TG2_scheme::loop_step_balance (tmesh::quadrant_iterator quadrant) //
{
  
  Dx = quadrant->p(0,1) - quadrant->p(0,0);
  Dy = quadrant->p(1,2) - quadrant->p(1,0);
  area = Dx*Dy;

  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){
      hdof[ii]    = sold_rkc [ordh    (quadrant->gt (ii))];
      Uxdof[ii]   = sold_rkc [ordUx   (quadrant->gt (ii))];
      Uydof[ii]   = sold_rkc [ordUy   (quadrant->gt (ii))];

      Z_node[ii] = Z[quadrant->gt (ii)];
      
      isdof_or_hanging[ii] = 1.;

    } else {
      hdof[ii]    = .5 * (sold_rkc [ordh  (quadrant->gparent(0,ii))] +
                          sold_rkc [ordh  (quadrant->gparent(1,ii))]);
      Uxdof[ii]   = .5 * (sold_rkc [ordUx (quadrant->gparent(0,ii))] +
                          sold_rkc [ordUx (quadrant->gparent(1,ii))]);
      Uydof[ii]   = .5 * (sold_rkc [ordUy (quadrant->gparent(0,ii))] +
                          sold_rkc [ordUy (quadrant->gparent(1,ii))]);

      Z_node[ii] = .5 * (Z [quadrant->gparent(0,ii)] +
                         Z [quadrant->gparent(1,ii)]);
      
      isdof_or_hanging[ii] = .5;

    }
  }

  const std::array<double,4> eta_vec = {hdof[0]+Z_node[0], hdof[1]+Z_node[1], hdof[2]+Z_node[2], hdof[3]+Z_node[3]};
  const std::array<double,4> der_x_eta = {(eta_vec[1]-eta_vec[0])/Dx, (eta_vec[1]-eta_vec[0])/Dx, (eta_vec[3]-eta_vec[2])/Dx, (eta_vec[3]-eta_vec[2])/Dx};
  const std::array<double,4> der_y_eta = {(eta_vec[2]-eta_vec[0])/Dy, (eta_vec[3]-eta_vec[1])/Dy, (eta_vec[2]-eta_vec[0])/Dy, (eta_vec[3]-eta_vec[1])/Dy};

  //std::cout << eta_vec[0] << " " << eta_vec[1] << " " << eta_vec[2] << " " << eta_vec[3] << std::endl;

  // trapezoidal rule integration!
  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){

      incr_source_balance [ordh  (quadrant->gt (ii))] += 0.;
      incr_source_balance [ordUx (quadrant->gt (ii))] += -grav*hdof[ii]*.25*area*der_x_eta[ii]*isdof_or_hanging[ii];
      incr_source_balance [ordUy (quadrant->gt (ii))] += -grav*hdof[ii]*.25*area*der_y_eta[ii]*isdof_or_hanging[ii];

    } else {

      incr_source_balance [ordh  (quadrant->gparent(0,ii))] += 0.; 
      incr_source_balance [ordh  (quadrant->gparent(1,ii))] += 0.;
      
      incr_source_balance [ordUx (quadrant->gparent(0,ii))] += -grav*hdof[ii]*.25*area*der_x_eta[ii]*isdof_or_hanging[ii];
      incr_source_balance [ordUx (quadrant->gparent(1,ii))] += -grav*hdof[ii]*.25*area*der_x_eta[ii]*isdof_or_hanging[ii];
      
      incr_source_balance [ordUy (quadrant->gparent(0,ii))] += -grav*hdof[ii]*.25*area*der_y_eta[ii]*isdof_or_hanging[ii];
      incr_source_balance [ordUy (quadrant->gparent(1,ii))] += -grav*hdof[ii]*.25*area*der_y_eta[ii]*isdof_or_hanging[ii];

    }
  }
  


}
*/



void
TG2_scheme::compute_stress_slope (tmesh::quadrant_iterator quadrant, const bool& isInitial)
{
  /*
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 
  const std::array<double,4> zeros_array = {0, 0, 0, 0};
  std::array<std::array<double,4>, 4 > delta_incr_mat = {zeros_array, zeros_array, zeros_array, zeros_array};

  std::array<int,4> bimpp_to_rev_ord = {0, 1, 3, 2};

  Dx = quadrant->p(0,1) - quadrant->p(0,0);
  Dy = quadrant->p(1,2) - quadrant->p(1,0);
  area = Dx*Dy;

  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){
      hdof[ii]    = sold_rkc_s [ordh    (quadrant->gt (ii))] + sold_rkc_w [ordh    (quadrant->gt (ii))];
      Uxdof[ii]   = sold_rkc_s [ordUx   (quadrant->gt (ii))] + sold_rkc_w [ordUx   (quadrant->gt (ii))];
      Uydof[ii]   = sold_rkc_s [ordUy   (quadrant->gt (ii))] + sold_rkc_w [ordUy   (quadrant->gt (ii))];
      
      isdof_or_hanging[ii] = 1.;

      delta_incr_mat[ii][ii] = tol_incr;
    } else {
      hdof[ii]    = .5 * (sold_rkc_s [ordh  (quadrant->gparent(0,ii))] +
                          sold_rkc_s [ordh  (quadrant->gparent(1,ii))]) +
                    .5 * (sold_rkc_w [ordh  (quadrant->gparent(0,ii))] +
                          sold_rkc_w [ordh  (quadrant->gparent(1,ii))]);
      Uxdof[ii]   = .5 * (sold_rkc_s [ordUx (quadrant->gparent(0,ii))] +
                          sold_rkc_s [ordUx (quadrant->gparent(1,ii))]) +
                    .5 * (sold_rkc_w [ordUx (quadrant->gparent(0,ii))] +
                          sold_rkc_w [ordUx (quadrant->gparent(1,ii))]);
      Uydof[ii]   = .5 * (sold_rkc_s [ordUy (quadrant->gparent(0,ii))] +
                          sold_rkc_s [ordUy (quadrant->gparent(1,ii))]) +
                    .5 * (sold_rkc_w [ordUy (quadrant->gparent(0,ii))] +
                          sold_rkc_w [ordUy (quadrant->gparent(1,ii))]);
      
      isdof_or_hanging[ii] = .5;

      delta_incr_mat[ii][ii] = tol_incr*.5;

      for (int jj = 0; jj < 4; ++jj)
      {
        if (quadrant->gparent(0,ii)==quadrant->gt (jj) || quadrant->gparent(1,ii)==quadrant->gt (jj))
        {
          delta_incr_mat[jj][ii] = tol_incr*.5;
        }

      }

    }
  }

  // weights coefficients for the flux term
  der_coeffs_x = {-Dy/2.*isdof_or_hanging[0], +Dy/2.*isdof_or_hanging[1],
    -Dy/2.*isdof_or_hanging[2], +Dy/2.*isdof_or_hanging[3]};
  
  der_coeffs_y = {-Dx/2.*isdof_or_hanging[0], -Dx/2.*isdof_or_hanging[1],
    +Dx/2.*isdof_or_hanging[2], +Dx/2.*isdof_or_hanging[3]};


  // compute the cell sigma_stress
  sigma_stress = is_stress_tensor ? compute_cell_stress (Uxdof[0], Uxdof[1], Uxdof[2], Uxdof[3], Uydof[0], Uydof[1], Uydof[2], Uydof[3]) : sigma_stress;

  for (int ii = 0; ii < 4; ++ii){
    D_U[ii] = U_stress_formula(hdof[ii], Uxdof[ii], Uydof[ii]);
  }

  std::array<std::array<double,3>,4> sigma_stress_incr_Ux, sigma_stress_incr_Uy;
  for (int ii = 0; ii < 4; ++ii){ // (3+3)x4

    sigma_stress_incr_Ux[ii] = is_stress_tensor ? compute_cell_stress (Uxdof[0]+delta_incr_mat[0][ii], Uxdof[1]+delta_incr_mat[1][ii], Uxdof[2]+delta_incr_mat[2][ii], Uxdof[3]+delta_incr_mat[3][ii], 
      Uydof[0], Uydof[1], Uydof[2], Uydof[3]) : std::array<double,3>{{0,0,0}};

    sigma_stress_incr_Uy[ii] = is_stress_tensor ? compute_cell_stress (Uxdof[0], Uxdof[1], Uxdof[2], Uxdof[3], 
      Uydof[0]+delta_incr_mat[0][ii], Uydof[1]+delta_incr_mat[1][ii], Uydof[2]+delta_incr_mat[2][ii], Uydof[3]+delta_incr_mat[3][ii]) : std::array<double,3>{{0,0,0}};
  }

  
  const std::array<double,4> der_phi_x = {-1./Dx*isdof_or_hanging[0], +1./Dx*isdof_or_hanging[1],
    -1./Dx*isdof_or_hanging[2], +1./Dx*isdof_or_hanging[3]};
  
  const std::array<double,4> der_phi_y = {-1./Dy*isdof_or_hanging[0], -1./Dy*isdof_or_hanging[1],
    +1./Dy*isdof_or_hanging[2], +1./Dy*isdof_or_hanging[3]};
  
    

  for (int ii = 0; ii < 4; ++ii){

    const double den1 = ii<2    ? 2. : 1.;
    const double den2 = ii<2    ? 1. : 2.;
    const double den3 = ii%2==1 ? 2. : 1.;
    const double den4 = ii%2==1 ? 1. : 2.;

    const double contribution_exact = (D_U[0]/den2+D_U[1]/den2+D_U[2]/den1+D_U[3]/den1);    

    
    const double h_  = 0.;
    const double Ux_ = area*.25*(der_phi_x[ii]*D_U[ii]*sigma_stress[0] + der_phi_y[ii]*D_U[ii]*sigma_stress[2]);
    const double Uy_ = area*.25*(der_phi_x[ii]*D_U[ii]*sigma_stress[2] + der_phi_y[ii]*D_U[ii]*sigma_stress[1]); 


    double h_s = 0., Ux_s = 0., Uy_s = 0.;
    for (int jj = 0; jj < 4; ++jj)
    {
      const auto sUx1_incr = (sigma_stress_incr_Ux[jj][0]-sigma_stress[0])/tol_incr;
      const auto sUx2_incr = (sigma_stress_incr_Ux[jj][2]-sigma_stress[2])/tol_incr;
      const auto sUx3_incr = (sigma_stress_incr_Uy[jj][0]-sigma_stress[0])/tol_incr;
      const auto sUx4_incr = (sigma_stress_incr_Uy[jj][2]-sigma_stress[2])/tol_incr;

      const auto & sUy1_incr = sUx2_incr;//(sigma_stress_incr_Ux[jj][2]-sigma_stress[2])/tol_incr;
      const auto   sUy2_incr = (sigma_stress_incr_Ux[jj][1]-sigma_stress[1])/tol_incr;
      const auto & sUy3_incr = sUx4_incr;//(sigma_stress_incr_Uy[jj][2]-sigma_stress[2])/tol_incr;
      const auto   sUy4_incr = (sigma_stress_incr_Uy[jj][1]-sigma_stress[1])/tol_incr;

      Ux_s += std::abs( area*.25*(der_phi_x[ii]*D_U[ii]*sUx1_incr + der_phi_y[ii]*D_U[ii]*sUx2_incr) );
      Ux_s += std::abs( area*.25*(der_phi_x[ii]*D_U[ii]*sUx3_incr + der_phi_y[ii]*D_U[ii]*sUx4_incr) );

      Uy_s += std::abs( area*.25*(der_phi_x[ii]*D_U[ii]*sUy1_incr + der_phi_y[ii]*D_U[ii]*sUy3_incr) );
      Uy_s += std::abs( area*.25*(der_phi_x[ii]*D_U[ii]*sUy3_incr + der_phi_y[ii]*D_U[ii]*sUy4_incr) );
    }
    

    if (isInitial)
    {
      if (! quadrant->is_hanging (ii)){

        stress_initial_step [ordh  (quadrant->gt (ii))] += h_;
        stress_initial_step [ordUx (quadrant->gt (ii))] += Ux_;
        stress_initial_step [ordUy (quadrant->gt (ii))] += Uy_;

        spec_radius_nodal [ordh   (quadrant->gt (ii))] += h_s;
        spec_radius_nodal [ordUx  (quadrant->gt (ii))] += Ux_s;
        spec_radius_nodal [ordUy  (quadrant->gt (ii))] += Uy_s;


      } else {

        stress_initial_step [ordh  (quadrant->gparent(0,ii))] += h_;
        stress_initial_step [ordh  (quadrant->gparent(1,ii))] += h_;

        stress_initial_step [ordUx (quadrant->gparent(0,ii))] += Ux_;
        stress_initial_step [ordUx (quadrant->gparent(1,ii))] += Ux_;

        stress_initial_step [ordUy (quadrant->gparent(0,ii))] += Uy_;
        stress_initial_step [ordUy (quadrant->gparent(1,ii))] += Uy_;



        spec_radius_nodal [ordh  (quadrant->gparent(0,ii))] += h_s;
        spec_radius_nodal [ordh  (quadrant->gparent(1,ii))] += h_s;

        spec_radius_nodal [ordUx (quadrant->gparent(0,ii))] += Ux_s;
        spec_radius_nodal [ordUx (quadrant->gparent(1,ii))] += Ux_s;

        spec_radius_nodal [ordUy (quadrant->gparent(0,ii))] += Uy_s;
        spec_radius_nodal [ordUy (quadrant->gparent(1,ii))] += Uy_s;

      }
    }
    else
    {
      if (! quadrant->is_hanging (ii)){

        stress_step [ordh  (quadrant->gt (ii))] += h_;
        stress_step [ordUx (quadrant->gt (ii))] += Ux_;
        stress_step [ordUy (quadrant->gt (ii))] += Uy_;

      } else {

        stress_step [ordh  (quadrant->gparent(0,ii))] += h_;
        stress_step [ordh  (quadrant->gparent(1,ii))] += h_;

        stress_step [ordUx (quadrant->gparent(0,ii))] += Ux_;
        stress_step [ordUx (quadrant->gparent(1,ii))] += Ux_;

        stress_step [ordUy (quadrant->gparent(0,ii))] += Uy_;
        stress_step [ordUy (quadrant->gparent(1,ii))] += Uy_;

      }
    }

  }
  */


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

    double hdof_c_s, Uxdof_c_s, Uydof_c_s, P_plus_h_c_s, P_minus_h_c_s, P_plus_Ux_c_s, P_minus_Ux_c_s, P_plus_Uy_c_s, P_minus_Uy_c_s;
    double hdof_c_w, Uxdof_c_w, Uydof_c_w, P_plus_h_c_w, P_minus_h_c_w, P_plus_Ux_c_w, P_minus_Ux_c_w, P_plus_Uy_c_w, P_minus_Uy_c_w;

    if (! quadrant->is_hanging (ii)){
      hdof_c_s      = sol_s [ordh    (quadrant->gt (ii))];
      Uxdof_c_s     = sol_s [ordUx   (quadrant->gt (ii))];
      Uydof_c_s     = sol_s [ordUy   (quadrant->gt (ii))];

      P_plus_h_c_s   = P_plus_s [ordh   (quadrant->gt (ii))];
      P_minus_h_c_s  = P_minus_s[ordh   (quadrant->gt (ii))];

      P_plus_Ux_c_s   = P_plus_s [ordUx   (quadrant->gt (ii))];
      P_minus_Ux_c_s  = P_minus_s[ordUx   (quadrant->gt (ii))];

      P_plus_Uy_c_s   = P_plus_s [ordUy   (quadrant->gt (ii))];
      P_minus_Uy_c_s  = P_minus_s[ordUy   (quadrant->gt (ii))];



      hdof_c_w      = sol_w [ordh    (quadrant->gt (ii))];
      Uxdof_c_w     = sol_w [ordUx   (quadrant->gt (ii))];
      Uydof_c_w     = sol_w [ordUy   (quadrant->gt (ii))];

      P_plus_h_c_w   = P_plus_w [ordh   (quadrant->gt (ii))];
      P_minus_h_c_w  = P_minus_w[ordh   (quadrant->gt (ii))];

      P_plus_Ux_c_w   = P_plus_w [ordUx   (quadrant->gt (ii))];
      P_minus_Ux_c_w  = P_minus_w[ordUx   (quadrant->gt (ii))];

      P_plus_Uy_c_w   = P_plus_w [ordUy   (quadrant->gt (ii))];
      P_minus_Uy_c_w  = P_minus_w[ordUy   (quadrant->gt (ii))];
      

    } else {
      hdof_c_s   = .5 * (sol_s [ordh  (quadrant->gparent(0,ii))] +
                         sol_s [ordh  (quadrant->gparent(1,ii))]);
      Uxdof_c_s  = .5 * (sol_s [ordUx (quadrant->gparent(0,ii))] +
                         sol_s [ordUx (quadrant->gparent(1,ii))]);
      Uydof_c_s  = .5 * (sol_s [ordUy (quadrant->gparent(0,ii))] +
                         sol_s [ordUy (quadrant->gparent(1,ii))]);

      P_plus_h_c_s   = .5 * (P_plus_s [ordh (quadrant->gparent(0,ii))] +
                             P_plus_s [ordh (quadrant->gparent(1,ii))]);
      P_minus_h_c_s  = .5 * (P_minus_s [ordh (quadrant->gparent(0,ii))] +
                             P_minus_s [ordh (quadrant->gparent(1,ii))]);

      P_plus_Ux_c_s   = .5 * (P_plus [ordUx (quadrant->gparent(0,ii))] +
                              P_plus [ordUx (quadrant->gparent(1,ii))]);
      P_minus_Ux_c_s  = .5 * (P_minus [ordUx (quadrant->gparent(0,ii))] +
                              P_mi_snus [ordUx (quadrant->gparent(1,ii))]);

      P_plus_Uy_c_s   = .5 * (P_plus_s [ordUy (quadrant->gparent(0,ii))] +
                              P_plus_s [ordUy (quadrant->gparent(1,ii))]);
      P_minus_Uy_c_s  = .5 * (P_minus_s [ordUy (quadrant->gparent(0,ii))] +
                              P_minus_s [ordUy (quadrant->gparent(1,ii))]);





      hdof_c_w   = .5 * (sol_w [ordh  (quadrant->gparent(0,ii))] +
                         sol_w [ordh  (quadrant->gparent(1,ii))]);
      Uxdof_c_w  = .5 * (sol_w [ordUx (quadrant->gparent(0,ii))] +
                         sol_w [ordUx (quadrant->gparent(1,ii))]);
      Uydof_c_w  = .5 * (sol_w [ordUy (quadrant->gparent(0,ii))] +
                         sol_w [ordUy (quadrant->gparent(1,ii))]);

      P_plus_h_c_w   = .5 * (P_plus_w [ordh (quadrant->gparent(0,ii))] +
                             P_plus_w [ordh (quadrant->gparent(1,ii))]);
      P_minus_h_c_w  = .5 * (P_minus_w [ordh (quadrant->gparent(0,ii))] +
                             P_minus_w [ordh (quadrant->gparent(1,ii))]);

      P_plus_Ux_c_w   = .5 * (P_plus_w [ordUx (quadrant->gparent(0,ii))] +
                              P_plus_w [ordUx (quadrant->gparent(1,ii))]);
      P_minus_Ux_c_w  = .5 * (P_minus_w [ordUx (quadrant->gparent(0,ii))] +
                              P_minus_w [ordUx (quadrant->gparent(1,ii))]);

      P_plus_Uy_c_w   = .5 * (P_plus_w [ordUy (quadrant->gparent(0,ii))] +
                              P_plus_w [ordUy (quadrant->gparent(1,ii))]);
      P_minus_Uy_c_w  = .5 * (P_minus_w [ordUy (quadrant->gparent(0,ii))] +
                              P_minus_w [ordUy (quadrant->gparent(1,ii))]);
      
    }

    hdof_s       [ii] = hdof_c_s;
    Uxdof_s      [ii] = Uxdof_c_s;
    Uydof_s      [ii] = Uydof_c_s;

    P_plus_h_dof_s [ii] = P_plus_h_c_s;
    P_minus_h_dof_s[ii] = P_minus_h_c_s;

    P_plus_Ux_dof_s [ii] = P_plus_Ux_c_s;
    P_minus_Ux_dof_s[ii] = P_minus_Ux_c_s;

    P_plus_Uy_dof_s [ii] = P_plus_Uy_c_s;
    P_minus_Uy_dof_s[ii] = P_minus_Uy_c_s;



    hdof_w       [ii] = hdof_c_w;
    Uxdof_w      [ii] = Uxdof_c_w;
    Uydof_w      [ii] = Uydof_c_w;

    P_plus_h_dof_w [ii] = P_plus_h_c_w;
    P_minus_h_dof_w[ii] = P_minus_h_c_w;

    P_plus_Ux_dof_w [ii] = P_plus_Ux_c_w;
    P_minus_Ux_dof_w[ii] = P_minus_Ux_c_w;

    P_plus_Uy_dof_w [ii] = P_plus_Uy_c_w;
    P_minus_Uy_dof_w[ii] = P_minus_Uy_c_w;


    hdof [ii]   = hdof_s [ii] + hdof_w [ii];
    Uxdof[ii]   = Uxdof_s[ii] + Uxdof_w[ii];
    Uydof[ii]   = Uydof_s[ii] + Uydof_w[ii];

  }
  
  
  
  // compute local extrema
  const auto h_min_cell_s  = *std::min_element(hdof_s.begin(),  hdof_s.end() );
  const auto h_max_cell_s  = *std::max_element(hdof_s.begin(),  hdof_s.end() );
 
  const auto Ux_min_cell_s = *std::min_element(Uxdof_s.begin(), Uxdof_s.end());
  const auto Ux_max_cell_s = *std::max_element(Uxdof_s.begin(), Uxdof_s.end());

  const auto Uy_min_cell_s = *std::min_element(Uydof_s.begin(), Uydof_s.end());
  const auto Uy_max_cell_s = *std::max_element(Uydof_s.begin(), Uydof_s.end()); 

  const auto h_min_cell_w  = *std::min_element(hdof_w.begin(),  hdof_w.end() );
  const auto h_max_cell_w  = *std::max_element(hdof_w.begin(),  hdof_w.end() );
 
  const auto Ux_min_cell_w = *std::min_element(Uxdof_w.begin(), Uxdof_w.end());
  const auto Ux_max_cell_w = *std::max_element(Uxdof_w.begin(), Uxdof_w.end());

  const auto Uy_min_cell_w = *std::min_element(Uydof_w.begin(), Uydof_w.end());
  const auto Uy_max_cell_w = *std::max_element(Uydof_w.begin(), Uydof_w.end());  

  bool is_node_in_element = false;

  std::array<double,4> h_min_s  = {h_min_cell_s, h_min_cell_s, h_min_cell_s, h_min_cell_s }, h_max_s  = {h_max_cell_s, h_max_cell_s, h_max_cell_s, h_max_cell_s },
                       Ux_min_s = {Ux_min_cell_s,Ux_min_cell_s,Ux_min_cell_s,Ux_min_cell_s}, Ux_max_s = {Ux_max_cell_s,Ux_max_cell_s,Ux_max_cell_s,Ux_max_cell_s},
                       Uy_min_s = {Uy_min_cell_s,Uy_min_cell_s,Uy_min_cell_s,Uy_min_cell_s}, Uy_max_s = {Uy_max_cell_s,Uy_max_cell_s,Uy_max_cell_s,Uy_max_cell_s};


  std::array<double,4> h_min_w  = {h_min_cell_w, h_min_cell_w, h_min_cell_w, h_min_cell_w }, h_max_w  = {h_max_cell_w, h_max_cell_w, h_max_cell_w, h_max_cell_w },
                       Ux_min_w = {Ux_min_cell_w,Ux_min_cell_w,Ux_min_cell_w,Ux_min_cell_w}, Ux_max_w = {Ux_max_cell_w,Ux_max_cell_w,Ux_max_cell_w,Ux_max_cell_w},
                       Uy_min_w = {Uy_min_cell_w,Uy_min_cell_w,Uy_min_cell_w,Uy_min_cell_w}, Uy_max_w = {Uy_max_cell_w,Uy_max_cell_w,Uy_max_cell_w,Uy_max_cell_w};



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

          double h_current_cell_s, Ux_current_cell_s, Uy_current_cell_s, h_current_cell_w, Ux_current_cell_w, Uy_current_cell_w;

          if (! quadrant_nei->is_hanging (jj)){
            h_current_cell_s  = sol_s [ordh  (quadrant_nei->gt (jj))];
            Ux_current_cell_s = sol_s [ordUx (quadrant_nei->gt (jj))];
            Uy_current_cell_s = sol_s [ordUy (quadrant_nei->gt (jj))];

            h_current_cell_w  = sol_w [ordh  (quadrant_nei->gt (jj))];
            Ux_current_cell_w = sol_w [ordUx (quadrant_nei->gt (jj))];
            Uy_current_cell_w = sol_w [ordUy (quadrant_nei->gt (jj))];
          } else {
            h_current_cell_s  = .5 * (sol_s [ordh  (quadrant_nei->gparent(0,jj))] +
                                      sol_s [ordh  (quadrant_nei->gparent(1,jj))]);
            Ux_current_cell_s = .5 * (sol_s [ordUx (quadrant_nei->gparent(0,jj))] +
                                      sol_s [ordUx (quadrant_nei->gparent(1,jj))]);
            Uy_current_cell_s = .5 * (sol_s [ordUy (quadrant_nei->gparent(0,jj))] +
                                      sol_s [ordUy (quadrant_nei->gparent(1,jj))]);

            h_current_cell_w  = .5 * (sol_w [ordh  (quadrant_nei->gparent(0,jj))] +
                                      sol_w [ordh  (quadrant_nei->gparent(1,jj))]);
            Ux_current_cell_w = .5 * (sol_w [ordUx (quadrant_nei->gparent(0,jj))] +
                                      sol_w [ordUx (quadrant_nei->gparent(1,jj))]);
            Uy_current_cell_w = .5 * (sol_w [ordUy (quadrant_nei->gparent(0,jj))] +
                                      sol_w [ordUy (quadrant_nei->gparent(1,jj))]);
          }

          h_min_s[ii]  = std::min(h_min_s[ii],  h_current_cell_s);
          h_max_s[ii]  = std::max(h_max_s[ii],  h_current_cell_s);

          Ux_min_s[ii] = std::min(Ux_min_s[ii], Ux_current_cell_s);
          Ux_max_s[ii] = std::max(Ux_max_s[ii], Ux_current_cell_s);

          Uy_min_s[ii] = std::min(Uy_min_s[ii], Uy_current_cell_s);
          Uy_max_s[ii] = std::max(Uy_max_s[ii], Uy_current_cell_s);


          h_min_w[ii]  = std::min(h_min_w[ii],  h_current_cell_w);
          h_max_w[ii]  = std::max(h_max_w[ii],  h_current_cell_w);

          Ux_min_w[ii] = std::min(Ux_min_w[ii], Ux_current_cell_w);
          Ux_max_w[ii] = std::max(Ux_max_w[ii], Ux_current_cell_w);

          Uy_min_w[ii] = std::min(Uy_min_w[ii], Uy_current_cell_w);
          Uy_max_w[ii] = std::max(Uy_max_w[ii], Uy_current_cell_w);

        }

        is_node_in_element = false;
      }


    }
  }



  
  // compute flux correction
  double phi_cell_h_s = 1., phi_cell_Ux_s = 1., phi_cell_Uy_s = 1.;
  double phi_cell_h_w = 1., phi_cell_Ux_w = 1., phi_cell_Uy_w = 1.;
  for (int ii = 0; ii < 4; ++ii){

    const auto & flux_on_the_node_h_s  = incr_anti_diff_s[ordh (index_quadrant)][ii];
    const auto & flux_on_the_node_Ux_s = incr_anti_diff_s[ordUx(index_quadrant)][ii];
    const auto & flux_on_the_node_Uy_s = incr_anti_diff_s[ordUy(index_quadrant)][ii];


    const auto & flux_on_the_node_h_w  = incr_anti_diff_w[ordh (index_quadrant)][ii];
    const auto & flux_on_the_node_Ux_w = incr_anti_diff_w[ordUx(index_quadrant)][ii];
    const auto & flux_on_the_node_Uy_w = incr_anti_diff_w[ordUy(index_quadrant)][ii];


    const auto & hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    const auto vel_rusanov_cell_x = hpoint>epsilon ? (std::abs(Uxdof[ii]/hpoint)+celerity) : 0.;
    const auto vel_rusanov_cell_y = hpoint>epsilon ? (std::abs(Uydof[ii]/hpoint)+celerity) : 0.;

    const auto vel_square_rusanov_cell = vel_rusanov_cell_x * vel_rusanov_cell_y;

    flux_limiter(h_min_s [ii], h_max_s [ii], hdof_s [ii], P_plus_h_dof_s [ii], P_minus_h_dof_s  [ii], flux_on_the_node_h_s,  vel_square_rusanov_cell, phi_cell_h_s );
    flux_limiter(Ux_min_s[ii], Ux_max_s[ii], Uxdof_s[ii], P_plus_Ux_dof_s[ii], P_minus_Ux_dof_s [ii], flux_on_the_node_Ux_s, vel_square_rusanov_cell, phi_cell_Ux_s);
    flux_limiter(Uy_min_s[ii], Uy_max_s[ii], Uydof_s[ii], P_plus_Uy_dof_s[ii], P_minus_Uy_dof_s [ii], flux_on_the_node_Uy_s, vel_square_rusanov_cell, phi_cell_Uy_s);


    flux_limiter(h_min_w [ii], h_max_w [ii], hdof_w [ii], P_plus_h_dof_w [ii], P_minus_h_dof_w  [ii], flux_on_the_node_h_w,  vel_square_rusanov_cell, phi_cell_h_w );
    flux_limiter(Ux_min_w[ii], Ux_max_w[ii], Uxdof_w[ii], P_plus_Ux_dof_w[ii], P_minus_Ux_dof_w [ii], flux_on_the_node_Ux_w, vel_square_rusanov_cell, phi_cell_Ux_w);
    flux_limiter(Uy_min_w[ii], Uy_max_w[ii], Uydof_w[ii], P_plus_Uy_dof_w[ii], P_minus_Uy_dof_w [ii], flux_on_the_node_Uy_w, vel_square_rusanov_cell, phi_cell_Uy_w);
  }


  for (int ii = 0; ii < 4; ++ii){

    const auto flux_on_the_node_h_s  = incr_anti_diff_s[ordh (index_quadrant)][ii]*phi_cell_h_s;
    const auto flux_on_the_node_Ux_s = incr_anti_diff_s[ordUx(index_quadrant)][ii]*phi_cell_Ux_s;
    const auto flux_on_the_node_Uy_s = incr_anti_diff_s[ordUy(index_quadrant)][ii]*phi_cell_Uy_s;


    const auto flux_on_the_node_h_w  = incr_anti_diff_w[ordh (index_quadrant)][ii]*phi_cell_h_w;
    const auto flux_on_the_node_Ux_w = incr_anti_diff_w[ordUx(index_quadrant)][ii]*phi_cell_Ux_w;
    const auto flux_on_the_node_Uy_w = incr_anti_diff_w[ordUy(index_quadrant)][ii]*phi_cell_Uy_w;

    if (! quadrant->is_hanging (ii)){

      incr_s [ordh  (quadrant->gt (ii))] += flux_on_the_node_h_s;
      incr_s [ordUx (quadrant->gt (ii))] += flux_on_the_node_Ux_s;
      incr_s [ordUy (quadrant->gt (ii))] += flux_on_the_node_Uy_s;

      incr_w [ordh  (quadrant->gt (ii))] += flux_on_the_node_h_w;
      incr_w [ordUx (quadrant->gt (ii))] += flux_on_the_node_Ux_w;
      incr_w [ordUy (quadrant->gt (ii))] += flux_on_the_node_Uy_w;

    } else {

      incr_s [ordh  (quadrant->gparent(0,ii))] += flux_on_the_node_h_s; 
      incr_s [ordh  (quadrant->gparent(1,ii))] += flux_on_the_node_h_s;
      
      incr_s [ordUx (quadrant->gparent(0,ii))] += flux_on_the_node_Ux_s;
      incr_s [ordUx (quadrant->gparent(1,ii))] += flux_on_the_node_Ux_s;
      
      incr_s [ordUy (quadrant->gparent(0,ii))] += flux_on_the_node_Uy_s;
      incr_s [ordUy (quadrant->gparent(1,ii))] += flux_on_the_node_Uy_s;


      incr_w [ordh  (quadrant->gparent(0,ii))] += flux_on_the_node_h_w; 
      incr_w [ordh  (quadrant->gparent(1,ii))] += flux_on_the_node_h_w;
      
      incr_w [ordUx (quadrant->gparent(0,ii))] += flux_on_the_node_Ux_w;
      incr_w [ordUx (quadrant->gparent(1,ii))] += flux_on_the_node_Ux_w;
      
      incr_w [ordUy (quadrant->gparent(0,ii))] += flux_on_the_node_Uy_w;
      incr_w [ordUy (quadrant->gparent(1,ii))] += flux_on_the_node_Uy_w;

    }

  }


}




void
TG2_scheme::flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& vel_square_rusanov_cell, double& phi_cell_Q)
{
  const auto Q_plus  = (Q_max-Q_dof)*(dt + dt_old)*.5*vel_square_rusanov_cell;
  const auto Q_minus = (Q_min-Q_dof)*(dt + dt_old)*.5*vel_square_rusanov_cell;

  const auto R_plus  = P_plus_Q ==0 ? 1 : std::min(1., Q_plus /P_plus_Q );
  const auto R_minus = P_minus_Q==0 ? 1 : std::min(1., Q_minus/P_minus_Q);

  phi_cell_Q = std::min(phi_cell_Q, flux_on_the_node>=0 ? R_plus : R_minus);

}

void
TG2_scheme::rkc(const int& j, const int& s, const int& kk)
{
  // kk is the current owned node

  double v_x_s, v_y_s, v_x_w, v_y_w;

  int count; 

  const double tolerance = 1.e-4;
  const int Nmax = 1e3;

  double error; 

  if (j == 1)
  {
    v_x_s = sol_ini_rkc_s.get_owned_data ()[kk+1] + mu_tilde_vect[1]*dt*stress_initial_step_s.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1];
    v_y_s = sol_ini_rkc_s.get_owned_data ()[kk+2] + mu_tilde_vect[1]*dt*stress_initial_step_s.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2];

    v_x_w = sol_ini_rkc_w.get_owned_data ()[kk+1] + mu_tilde_vect[1]*dt*stress_initial_step_w.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1];
    v_y_w = sol_ini_rkc_w.get_owned_data ()[kk+2] + mu_tilde_vect[1]*dt*stress_initial_step_w.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2];
  }
  else
  {
    v_x_s = (1. - mu_vect[j] - v_vect[j])*sol_ini_rkc_s.get_owned_data ()[kk+1] + mu_vect[j]*sold_rkc_s.get_owned_data ()[kk+1] + 
    v_vect[j]*soldd_rkc_s.get_owned_data ()[kk+1] + mu_tilde_vect[j]*dt*stress_step_s.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1] - v_vect[j]*mu_tilde_vect[1]*dt*incr_source_s.get_owned_data ()[kk+1];
    gamma_tilde_vect[j]*dt*stress_initial_step_s.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1] + 
    (gamma_tilde_vect[j]*mu_tilde_vect[1]*mu_vect[j]/mu_tilde_vect[j] - (1. - mu_vect[j] - v_vect[j])*mu_tilde_vect[1])*dt*incr_initial_source_s.get_owned_data ()[kk+1] - v_vect[j]*mu_tilde_vect[1]*dt*incr_source_s.get_owned_data ()[kk+1];

    v_y_s = (1. - mu_vect[j] - v_vect[j])*sol_ini_rkc_s.get_owned_data ()[kk+2] + mu_vect[j]*sold_rkc_s.get_owned_data ()[kk+2] + 
    v_vect[j]*soldd_rkc_s.get_owned_data ()[kk+2] + mu_tilde_vect[j]*dt*stress_step_s.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] + 
    gamma_tilde_vect[j]*dt*stress_initial_step_s.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] +
    (gamma_tilde_vect[j]*mu_tilde_vect[1]*mu_vect[j]/mu_tilde_vect[j] - (1. - mu_vect[j] - v_vect[j])*mu_tilde_vect[1])*dt*incr_initial_source_s.get_owned_data ()[kk+2] - v_vect[j]*mu_tilde_vect[1]*dt*incr_source_s.get_owned_data ()[kk+2];



    v_x_w = (1. - mu_vect[j] - v_vect[j])*sol_ini_rkc_w.get_owned_data ()[kk+1] + mu_vect[j]*sold_rkc_w.get_owned_data ()[kk+1] + 
    v_vect[j]*soldd_rkc_w.get_owned_data ()[kk+1] + mu_tilde_vect[j]*dt*stress_step_w.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1] - v_vect[j]*mu_tilde_vect[1]*dt*incr_source_w.get_owned_data ()[kk+1];
    gamma_tilde_vect[j]*dt*stress_initial_step_w.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1] + 
    (gamma_tilde_vect[j]*mu_tilde_vect[1]*mu_vect[j]/mu_tilde_vect[j] - (1. - mu_vect[j] - v_vect[j])*mu_tilde_vect[1])*dt*incr_initial_source_w.get_owned_data ()[kk+1] - v_vect[j]*mu_tilde_vect[1]*dt*incr_source_w.get_owned_data ()[kk+1];

    v_y_w = (1. - mu_vect[j] - v_vect[j])*sol_ini_rkc_w.get_owned_data ()[kk+2] + mu_vect[j]*sold_rkc_w.get_owned_data ()[kk+2] + 
    v_vect[j]*soldd_rkc_w.get_owned_data ()[kk+2] + mu_tilde_vect[j]*dt*stress_step_w.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] + 
    gamma_tilde_vect[j]*dt*stress_initial_step_w.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] +
    (gamma_tilde_vect[j]*mu_tilde_vect[1]*mu_vect[j]/mu_tilde_vect[j] - (1. - mu_vect[j] - v_vect[j])*mu_tilde_vect[1])*dt*incr_initial_source_w.get_owned_data ()[kk+2] - v_vect[j]*mu_tilde_vect[1]*dt*incr_source_w.get_owned_data ()[kk+2];
  }


  // solve non-linearities
  count = -1;
  error = tolerance + 1; 

  const auto & h_c_s  = sol_s.   get_owned_data ()[    kk  ];
  const auto & h_c_w  = sol_w.   get_owned_data ()[    kk  ];
  const auto & poro_c = porosity.get_owned_data ()[int(kk/3)];

  const auto h_c = h_c_w + h_c_s;

  const double C_w = grav*(density_s-density_w)/terminal_velocity/std::pow(poro_c,m_coeff)/density_w;
  const double C_s = grav*(density_s-density_w)/terminal_velocity/std::pow(poro_c,m_coeff)/density_s;

  const double densi_ = (1-poro_c)*(density_s - density_w)
  while (count++<Nmax && error>tolerance)
  {
    const auto & Ux_c_s = sol_s.get_owned_data ()[kk+1];
    const auto & Uy_c_s = sol_s.get_owned_data ()[kk+2];

    const auto & Ux_c_w = sol_w.get_owned_data ()[kk+1];
    const auto & Uy_c_w = sol_w.get_owned_data ()[kk+2];

    const auto U_tot_x = Ux_c_s + Ux_c_w;
    const auto U_tot_y = Uy_c_s + Uy_c_w;

    const auto abs_mass_flux = std::sqrt(U_tot_x*U_tot_x + U_tot_y*U_tot_y);


    const auto big_A   = (1. - mu_tilde_vect[1]*dt*C_w*(1.-poro_c));
    const auto big_B   = mu_tilde_vect[1]*dt*C_w*poro_c;
    const auto big_C_x = - mu_tilde_vect[1]*dt*C_s*(1.-poro_c) - mu_tilde_vect[1]*dt/density_s*( (h_c*h_c>epsilon && is_bed_friction) ? -density*grav/h_c/h_c/turbulence_coeff*2.*std::abs(U_tot_x) : 0. + is_bed_friction ? -densi_*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign) : 0. );
    const auto big_C_y = - mu_tilde_vect[1]*dt*C_s*(1.-poro_c) - mu_tilde_vect[1]*dt/density_s*( (h_c*h_c>epsilon && is_bed_friction) ? -density*grav/h_c/h_c/turbulence_coeff*2.*std::abs(U_tot_y) : 0. + is_bed_friction ? -densi_*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign) : 0. );
    const auto big_D_x = 1. + mu_tilde_vect[1]*dt*C_s*poro_c   - mu_tilde_vect[1]*dt/density_s*( (h_c*h_c>epsilon && is_bed_friction) ? -density*grav/h_c/h_c/turbulence_coeff*2.*std::abs(U_tot_x) : 0. + is_bed_friction ? -densi_*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign) : 0. );
    const auto big_D_y = 1. + mu_tilde_vect[1]*dt*C_s*poro_c   - mu_tilde_vect[1]*dt/density_s*( (h_c*h_c>epsilon && is_bed_friction) ? -density*grav/h_c/h_c/turbulence_coeff*2.*std::abs(U_tot_y) : 0. + is_bed_friction ? -densi_*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign) : 0. );


    const auto F_Ux_w = v_x_w - Ux_c_w + mu_tilde_vect[1]*dt*Ux_src_formula_w(h_c, Ux_c_w, Ux_c_s, Uy_c_w, Uy_c_s, poro_c);
    const auto F_Ux_s = v_x_s - Ux_c_s + mu_tilde_vect[1]*dt*Ux_src_formula_s(h_c, Ux_c_w, Ux_c_s, Uy_c_w, Uy_c_s, poro_c);
    const auto F_Uy_w = v_y_w - Uy_c_w + mu_tilde_vect[1]*dt*Uy_src_formula_w(h_c, Ux_c_w, Ux_c_s, Uy_c_w, Uy_c_s, poro_c);
    const auto F_Uy_s = v_y_s - Uy_c_s + mu_tilde_vect[1]*dt*Uy_src_formula_s(h_c, Ux_c_w, Ux_c_s, Uy_c_w, Uy_c_s, poro_c);



    const auto delta_Ux_s = (F_Ux_s - big_C_x/big_A*F_Ux_w)/(big_D_x - big_C_x/big_A*big_B);
    const auto delta_Uy_s = (F_Uy_s - big_C_y/big_A*F_Uy_w)/(big_D_y - big_C_y/big_A*big_B);

    const auto delta_Ux_w = (F_Ux_w - big_B*delta_Ux_s)/big_A;
    const auto delta_Uy_w = (F_Uy_w - big_B*delta_Uy_s)/big_A;

    error = std::sqrt(delta_Ux_s*delta_Ux_s + delta_Uy_s*delta_Uy_s + delta_Ux_w*delta_Ux_w + delta_Uy_w*delta_Uy_w);


    sol_s.get_owned_data ()[kk+1] += delta_Ux_s; 
    sol_w.get_owned_data ()[kk+1] += delta_Ux_w;

    sol_s.get_owned_data ()[kk+2] += delta_Uy_s; 
    sol_w.get_owned_data ()[kk+2] += delta_Uy_w;

  } 
  if (j == s) Newton_it.get_owned_data ()[int(kk/3)] = double(count);



}


void
TG2_scheme::prepare_IMEXRKC_coefficients (const int& s)
{
  // \omega_0
  w0 = 1. + epsilon_IMEXRKC/s/s;

  double T, T_old = w0, T_oldold = 1., 
  T_prime, T_prime_old = 1., T_prime_oldold = 0.,
  T_second, T_second_old = 0., T_second_oldold = 0.; 


  b_vect.resize(s+1);
  gamma_tilde_vect.resize(s+1);
  for (int iii=2; iii<=s; iii++)
  {
    // first Chebyshev function
    T = 2.*w0*T_old - T_oldold;

    // first Chebyshev function derivative
    T_prime = 2.*T_old + 2.*w0*T_prime_old - T_prime_oldold;

    // first Chebyshev function second derivative 
    T_second = 4.*T_prime_old + 2.*w0*T_second_old - T_second_oldold;


    b_vect[iii] = T_second/T_prime/T_prime;

    // -a_{j-1} contribution
    gamma_tilde_vect[iii] = -(1. - b_vect[iii-1]*T_old);

    // update
    T_oldold = T_old;
    T_old = T;

    T_prime_oldold = T_prime_old;
    T_prime_old = T_prime;

    T_second_oldold = T_second_old;
    T_second_old = T_second;
  }

  b_vect[0] = b_vect[2];
  b_vect[1] = b_vect[2];


  // -a_{0} contribution
  gamma_tilde_vect[1] = -(1. - b_vect[0]*w0);

  // \omega_1
  w1 = T_prime/T_second;

  mu_tilde_vect.resize(s+1);
  v_vect.resize(s+1);
  mu_vect.resize(s+1);


  mu_tilde_vect[0] = 0.;
  mu_tilde_vect[1] = (s == 1) ? 1. : b_vect[1]*w1;

  gamma_tilde_vect[1] = gamma_tilde_vect[1]*mu_tilde_vect[1];
  for (int iii=2; iii<=s; iii++)
  {
    mu_tilde_vect[iii] = 2.*b_vect[iii]*w1/b_vect[iii-1];

    gamma_tilde_vect[iii] = gamma_tilde_vect[iii]*mu_tilde_vect[iii];

    v_vect[iii] = -b_vect[iii]/b_vect[iii-2];

    mu_vect[iii] = 2.*b_vect[iii]*w0/b_vect[iii-1];

    //std::cout << mu_tilde_vect[iii] << " " << w1  << " " << gamma_tilde_vect[iii] << " " << v_vect[iii] << " " << mu_vect[iii] << " " << b_vect[iii] << std::endl;
  
  }


/*
  w_fun_0(s);

  T_fun       (s);
  T_fun_prime (s);
  T_fun_second(s);

  w_fun_1(s);

  // I want just these to take in memory, 5 arrays of storage!!
  b_fun(s);
  mu_fun_tilde(s);
  gamma_tilde_fun(s);
  v_fun(s);
  mu_fun(s);*/



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
{ return Ux; }

double
TG2_scheme::h_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return Uy; }

double
TG2_scheme::Ux_flux_formula_x (const double& h, const double& Ux, const double& Uy, const double& h_)
{ 
  const auto vel_x = h>epsilon ? Ux/h : 0.;
  return (Ux*vel_x + grav*h*(h+h_)/2.); 
}
 
double
TG2_scheme::Ux_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return (h>epsilon ? Uy*Ux/h : 0.); }

double
TG2_scheme::Uy_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ return (h>epsilon ? Uy*Ux/h : 0.); }

double
TG2_scheme::Uy_flux_formula_y (const double& h, const double& Ux, const double& Uy, const double& h_)
{ 
  const auto vel_y = h>epsilon ? Uy/h : 0.;
  return (Uy*vel_y + grav*h*(h+h_)/2.); 
}


// stress functions
double
TG2_scheme::U_stress_formula (const double& h, const double& Ux, const double& Uy)
{ return (-h/density); }




std::array<double,3>
TG2_scheme::compute_cell_stress (const double& Uxdof_0, const double& Uxdof_1, 
                                 const double& Uxdof_2, const double& Uxdof_3, 
                                 const double& Uydof_0, const double& Uydof_1, 
                                 const double& Uydof_2, const double& Uydof_3)
{
  // compute \sigma_xx, ...

  // def_grad = [D11, D22, D33, D12, D23, D31]
  // sigma = [sigma_11, sigma_22, sigma_12]

  //std::cout << Uxdof_0 << " " << Uxdof_1 << " " << Uxdof_2 << " " << Uxdof_3 << " " << Uydof_0 << " " << Uydof_1 << " " << Uydof_2 << " " << Uydof_3 << std::endl;

  std::array<double,6> def_grad = compute_cell_def_grad (Uxdof_0, Uxdof_1, Uxdof_2, Uxdof_3, 
                                                         Uydof_0, Uydof_1, Uydof_2, Uydof_3);

  double second_invariant = 0.;
  for (int i_def = 0; i_def < 6; i_def++)
  {
    second_invariant += def_grad[i_def]*def_grad[i_def];
  }
  second_invariant *= .5;

  const double kinetic_energergy_associated = std::sqrt(second_invariant);

  double viscos = kinetic_energergy_associated!=0 ? 2.*fluid_viscosity + yield_shear_stress/kinetic_energergy_associated*(1. - std::exp(-regularization_parameter*kinetic_energergy_associated)) : 0.;
  viscos = yield_shear_stress==0 ? 2.*fluid_viscosity : viscos;

  //std::cout << def_grad[0] << " " << def_grad[1] << " " << def_grad[2] << " " << def_grad[3] << " " << def_grad[4] << " " << def_grad[5] << " aa" << std::endl;

//std::cout << viscos << std::endl;//def_grad[0] << " " << def_grad[1] << " " << def_grad[2] << std::endl;
  return(std::array<double,3>{{viscos*def_grad[0], viscos*def_grad[1], viscos*def_grad[3]}});
}


std::array<double,6>
TG2_scheme::compute_cell_def_grad (const double& Uxdof_0, const double& Uxdof_1, 
                                   const double& Uxdof_2, const double& Uxdof_3, 
                                   const double& Uydof_0, const double& Uydof_1, 
                                   const double& Uydof_2, const double& Uydof_3)
{
  // def_grad = [D11, D22, D33, D12, D23, D31]

  // compute \zeta
  double h_cell = .25*(hdof[0]+hdof[1]+hdof[2]+hdof[3]), 
  Ux_cell = .25*(Uxdof_0+Uxdof_1+Uxdof_2+Uxdof_3), Uy_cell = .25*(Uydof_0+Uydof_1+Uydof_2+Uydof_3); 

  //std::cout << Uxdof_0 << " " << hdof[0] << std::endl;

  const auto & vel_x_3 = hdof[3]>epsilon ? Uxdof_3/hdof[3] : 0.;
  const auto & vel_x_2 = hdof[2]>epsilon ? Uxdof_2/hdof[2] : 0.;
  const auto & vel_x_1 = hdof[1]>epsilon ? Uxdof_1/hdof[1] : 0.;
  const auto & vel_x_0 = hdof[0]>epsilon ? Uxdof_0/hdof[0] : 0.;

  const auto & vel_y_3 = hdof[3]>epsilon ? Uydof_3/hdof[3] : 0.;
  const auto & vel_y_2 = hdof[2]>epsilon ? Uydof_2/hdof[2] : 0.;
  const auto & vel_y_1 = hdof[1]>epsilon ? Uydof_1/hdof[1] : 0.;
  const auto & vel_y_0 = hdof[0]>epsilon ? Uydof_0/hdof[0] : 0.;

  grad_cell_ux = {.5 * ( (vel_x_3 - vel_x_2) + (vel_x_1 - vel_x_0) )/Dx, .5 * ( (vel_x_2 - vel_x_0) + (vel_x_3 - vel_x_1) )/Dy};
  grad_cell_uy = {.5 * ( (vel_y_3 - vel_y_2) + (vel_y_1 - vel_y_0) )/Dx, .5 * ( (vel_y_2 - vel_y_0) + (vel_y_3 - vel_y_1) )/Dy};



  const double vel_x_cell = h_cell>epsilon ? Ux_cell/h_cell : 0.;
  const double vel_y_cell = h_cell>epsilon ? Uy_cell/h_cell : 0.;


  const double abs_vel_cell = std::sqrt( vel_x_cell*vel_x_cell + vel_y_cell*vel_y_cell );

  const double a = 3./2.;
  const double c = 65./32.;
  const double aa = h_cell>epsilon ? 6.*fluid_viscosity*abs_vel_cell/h_cell/yield_shear_stress : 0.;


  const double b = -(114./32.+aa);
  const double Delta = b*b - 4.*a*c;

  const double zeta_1 = (-b + std::sqrt(Delta))/2./a;
  const double zeta_2 = (-b - std::sqrt(Delta))/2./a;

  if ( std::abs(zeta_1 - .5)<=.5 && std::abs(zeta_2 - .5)<=.5)
  {
    std::cout << "Two valid roots, look at compute_cell_def_grad function, " << zeta_1 << " " << zeta_2 << ", STOP!" << std::endl;
    exit(1.);
  }

  const double zeta = std::abs(zeta_1 - .5)<=.5 ? zeta_1 : zeta_2;

  //std::cout << zeta << " " << aa << " " << vel_x_cell << " " << vel_y_cell << " " << abs_vel_cell << " bbbb" << std::endl;


  const auto & partial_x_ux = grad_cell_ux[0];
  const auto & partial_y_ux = grad_cell_ux[1];
  const auto   partial_z_ux = h_cell>epsilon ? 3./(2.+zeta)*vel_x_cell/h_cell : 0.;

  const auto & partial_x_uy = grad_cell_uy[0];
  const auto & partial_y_uy = grad_cell_uy[1];
  const auto   partial_z_uy = h_cell>epsilon ? 3./(2.+zeta)*vel_y_cell/h_cell : 0.;

  const auto partial_x_uz = 0.; // steady state simple shear flow
  const auto partial_y_uz = 0.; // steady state simple shear flow
  const auto partial_z_uz = -(grad_cell_ux[0]+grad_cell_uy[1]);
  
  return(std::array<double,6>{{partial_x_ux, partial_y_uy, partial_z_uz, .5*(partial_x_uy+partial_y_ux), .5*(partial_z_uy+partial_y_uz), .5*(partial_z_ux+partial_x_uz)}});
}


// source terms
double
TG2_scheme::h_src_formula (const double& h, const double& Ux, const double& Uy, const double& n_poro)
{ 
  const double abs_U = std::sqrt( Ux*Ux + Uy*Uy );
  return (is_erosion ? n_poro*erosion_coefficient*abs_U : 0.); 
}

double
TG2_scheme::Ux_src_formula_w (const double& h, const double& U_wx, const double& U_sx, const double& U_wy, const double& U_sy, const double& n_poro)
{ 
  const auto contr_2 = int_term_x(h, U_wx, U_sx, n_poro, density_w);

  return(contr_2);
}

double
TG2_scheme::Ux_src_formula_s (const double& h, const double& U_wx, const double& U_sx, const double& U_wy, const double& U_sy, const double& n_poro)
{ 
  const auto contr_1 = friction_Ux (h, Ux_wx+U_sx, U_wy+U_sy, n_poro);
  const auto contr_2 = int_term_x(h, U_wx, U_sx, n_poro, density_s);

  return(contr_1 + contr_2);
}

double
TG2_scheme::Uy_src_formula_w (const double& h, const double& U_wx, const double& U_sx, const double& U_wy, const double& U_sy, const double& n_poro)
{ 
  const auto contr_2 = int_term_y(h, U_wy, U_sy, n_poro, density_w);

  return(contr_2);
}

double
TG2_scheme::Uy_src_formula_s (const double& h, const double& U_wx, const double& U_sx, const double& U_wy, const double& U_sy, const double& n_poro)
{ 
  const auto contr_1 = friction_Uy (h, Ux_wx+U_sx, U_wy+U_sy, n_poro);
  const auto contr_2 = int_term_y(h, U_wy, U_sy, n_poro, density_s);

  return(contr_1 + contr_2);
}

double
TG2_scheme::int_term_y (const double& h, const double& U_wy, const double& U_sy, const double& n_poro, const double& dens_a)
{
  const double C_a = grav/dens_a*(density_s-density_w) / (terminal_velocity*std::pow(n_poro,m_coeff));
  return( C_a*((1.-n_poro)*U_wy - n_poro*U_sy) );
}

double
TG2_scheme::int_term_x (const double& h, const double& U_wx, const double& U_sx, const double& n_poro, const double& dens_a)
{
  const double C_a = grav/dens_a*(density_s-density_w) / (terminal_velocity*std::pow(n_poro,m_coeff));
  return( C_a*((1.-n_poro)*U_wx - n_poro*U_sx) );
}

double
TG2_scheme::friction_Ux (const double& h, const double& Ux, const double& Uy, const double& n_poro)
{
  const double densi_ = (1-n_poro)*(density_s - density_w);
  const double bed_pressure = densi_*grav*h/density_s; //grav*h + surface_pressure/density; 

  const double vel_x_sign = std::abs(Ux)>tolerance_sign ? Ux/std::abs(Ux) : Ux/tolerance_sign;

  const double bed_fric_contr_one = (h*h>epsilon && is_bed_friction) ? Ux*grav*std::abs(Ux)/turbulence_coeff/h/h/density_s*density : 0.; 
  const double bed_fric_contr_two = is_bed_friction ? vel_x_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;

  return (- bed_fric_contr_one - bed_fric_contr_two);
}

double
TG2_scheme::friction_Uy (const double& h, const double& Ux, const double& Uy, const double& n_poro)
{
  const double densi_ = (1-n_poro)*(density_s - density_w);
  const double bed_pressure = densi_*grav*h/density_s; // grav*h + surface_pressure/density;

  const double vel_y_sign = std::abs(Uy)>tolerance_sign ? Uy/std::abs(Uy) : Uy/tolerance_sign;

  const double bed_fric_contr_one = (h*h>epsilon && is_bed_friction) ? Uy*grav*std::abs(Uy)/turbulence_coeff/h/h/density_s*density  : 0.; //is_bed_friction ? vel_y*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = is_bed_friction ? vel_y_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;

  return (- bed_fric_contr_one - bed_fric_contr_two);
}


double
TG2_scheme::src_slope_formula (const double& h, const double& S_x, const double& S_y, const int& kk)
{
  // nodal slope term for the second step in the main

  double S = 0.;

  if (kk%3==1) // Ux
  {
    S = S_x;
  }
  else if (kk%3==2) // Uy
  {
    S = S_y;
  }

  //std::cout << S << " " << kk << " " << kk%3 << std::endl;

  return (-grav*S*h);

}


double 
TG2_scheme::src_slope_formula (const double& h, const double& S)
{
  // cell-wise source term to build the incr vector
  return (-grav*S*h);
}





