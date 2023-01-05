#include "Taylor_Galerkin_IMEX-RKC_Strang_balanced_two_phase.h"
#include <algorithm>
#include <cassert>

TG2_scheme::TG2_scheme(Q1& sol,
                       Q1& sold,
                       Q1& soldd, 
                       Q1& sold_rkc,
                       Q1& soldd_rkc, 
                       Q1& sol_ini_rkc,
                       Q1& incr,
                       Q1& incr_initial_source,
                       Q1& incr_source,
                       std::vector<std::array<double,4>>& incr_anti_diff,
                       Q1& stress_initial_step,
                       Q1& stress_step,
                       Q1& P_plus,
                       Q1& P_minus,
                       Q1& spec_radius_nodal,
                       Q0& sol_onehalf,
                       Q1& mass,
                       const ordering& oh,
                       const ordering& on,
                       const ordering& oUxw,
                       const ordering& oUyw,
                       const ordering& oUxs,
                       const ordering& oUys,
                       const Q1& Z,
                       Q0& Z_onehalf,
		                   Q1& Newton_it,
                       const double& DELTAT,
                       const double& h_min,
                       const bool& is_non_reflBC,
                       const bool& is_stress_tensor,
                       const double& grav,
                       const double& density_w,
                       const double& density_s,
                       const double& turbulence_coeff,
                       const double& bed_friction_angle_rad,
                       const double& fluid_viscosity,
                       const double& yield_shear_stress,
                       const double& erosion_coefficient,
                       const double& m_coeff,
                       const double& terminal_velocity)
: sol(sol), sold(sold), soldd(soldd), sold_rkc(sold_rkc), soldd_rkc(soldd_rkc), sol_ini_rkc(sol_ini_rkc), incr(incr), incr_initial_source(incr_initial_source), incr_source(incr_source), incr_anti_diff(incr_anti_diff), stress_initial_step(stress_initial_step), stress_step(stress_step), P_plus(P_plus), P_minus(P_minus), spec_radius_nodal(spec_radius_nodal), sol_onehalf(sol_onehalf), mass(mass), 
  ordh(oh), ordn(on), ordUxw(oUxw), ordUyw(oUyw), ordUxs(oUxs), ordUys(oUys), Z(Z), Z_onehalf(Z_onehalf), Newton_it(Newton_it), DELTAT(DELTAT), epsilon(h_min), is_non_reflBC(is_non_reflBC), is_stress_tensor(is_stress_tensor), grav(grav),
  density_w(density_w), density_s(density_s), turbulence_coeff(turbulence_coeff), bed_friction_angle_rad(bed_friction_angle_rad), fluid_viscosity(fluid_viscosity), yield_shear_stress(yield_shear_stress), erosion_coefficient(erosion_coefficient), 
  m_coeff(m_coeff), terminal_velocity(terminal_velocity)
{ }
 

std::array<double,2>
TG2_scheme::max_eigen (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  double lambda_x = 0., 
         lambda_y = 0., 
         hw = n*h,
         hs = (1.-n)*h,

         vel_x  = h>epsilon ? (Uxw+Uxs)/h : 0., 
         vel_y  = h>epsilon ? (Uyw+Uys)/h : 0., 
         velw_x = hw>epsilon ? Uxw/hw : 0., 
         velw_y = hw>epsilon ? Uyw/hw : 0., 
         vels_x = hs>epsilon ? Uxs/hs : 0., 
         vels_y = hs>epsilon ? Uys/hs : 0., 

         celerity = std::sqrt(grav*h),
         g_prime = (1.-density_w/density_s)*grav;


  lambda_x = h>epsilon ? std::max(std::abs(vel_x)+celerity, std::abs(velw_x*hs+vels_x*hw/h)+std::sqrt(hw*hs/h*(g_prime-(velw_x-vels_x)*(velw_x-vels_x)/h) )) : 0.;
  lambda_y = h>epsilon ? std::max(std::abs(vel_y)+celerity, std::abs(velw_y*hs+vels_y*hw/h)+std::sqrt(hw*hs/h*(g_prime-(velw_y-vels_y)*(velw_y-vels_y)/h) )) : 0.;

  return(std::array<double,2>{{lambda_x, lambda_y}});
}

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
      hdof[ii]   = sol [ordh   (quadrant->gt (ii) )];
      ndof[ii]   = sol [ordn   (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 
    }
    else
    {
      hdof[ii]   = .5 * (sol [ordh   (quadrant->gparent (0, ii) )] +
                         sol [ordh   (quadrant->gparent (1, ii) )]);
      ndof[ii]   = .5 * (sol [ordn   (quadrant->gparent (0, ii) )] +
                         sol [ordn   (quadrant->gparent (1, ii) )]);
      Uxwdof[ii] = .5 * (sol [ordUxw (quadrant->gparent (0, ii) )] +
                         sol [ordUxw (quadrant->gparent (1, ii) )]);
      Uywdof[ii] = .5 * (sol [ordUyw (quadrant->gparent (0, ii) )] +
                         sol [ordUyw (quadrant->gparent (1, ii) )]);
      Uxsdof[ii] = .5 * (sol [ordUxs (quadrant->gparent (0, ii) )] +
                         sol [ordUxs (quadrant->gparent (1, ii) )]);
      Uysdof[ii] = .5 * (sol [ordUys (quadrant->gparent (0, ii) )] +
                         sol [ordUys (quadrant->gparent (1, ii) )]);
    }
  }

  
  for (int ii = 0; ii < 4; ++ii){

    const auto lambdas = max_eigen (hdof[ii], ndof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii]);
    const auto & vel_rusanov_cell_x = lambdas[0];
    const auto & vel_rusanov_cell_y = lambdas[1];

    const auto dtoptx = hdof[ii]>epsilon ? Dx/vel_rusanov_cell_x : DELTAT;
    const auto dtopty = hdof[ii]>epsilon ? Dy/vel_rusanov_cell_y : DELTAT;
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
  
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 
  
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii);
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;


  double h_cell_average = 0., n_cell_average = 0., Uwx_cell_average = 0., Uwy_cell_average = 0., Usx_cell_average = 0., Usy_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hdof[ii]   = sol [ordh   (quadrant->gt (ii) )];
      ndof[ii]   = sol [ordn   (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 

      Z_node[ii]  = Z [quadrant->gt (ii)];
    }
    else
    {
      hdof[ii]   = .5 * (sol [ordh   (quadrant->gparent (0, ii) )] +
                         sol [ordh   (quadrant->gparent (1, ii) )]);
      ndof[ii]   = .5 * (sol [ordn   (quadrant->gparent (0, ii) )] +
                         sol [ordn   (quadrant->gparent (1, ii) )]);
      Uxwdof[ii] = .5 * (sol [ordUxw (quadrant->gparent (0, ii) )] +
                         sol [ordUxw (quadrant->gparent (1, ii) )]);
      Uywdof[ii] = .5 * (sol [ordUyw (quadrant->gparent (0, ii) )] +
                         sol [ordUyw (quadrant->gparent (1, ii) )]);
      Uxsdof[ii] = .5 * (sol [ordUxs (quadrant->gparent (0, ii) )] +
                         sol [ordUxs (quadrant->gparent (1, ii) )]);
      Uysdof[ii] = .5 * (sol [ordUys (quadrant->gparent (0, ii) )] +
                         sol [ordUys (quadrant->gparent (1, ii) )]);

      Z_node[ii]  = .5 * (Z [quadrant->gparent(0,ii)] +
                          Z [quadrant->gparent(1,ii)]);
    }

    const auto & hdof_c   = hdof  [ii];
    const auto & ndof_c   = ndof  [ii];
    const auto & Uxwdof_c = Uxwdof[ii];
    const auto & Uywdof_c = Uywdof[ii];
    const auto & Uxsdof_c = Uxsdof[ii];
    const auto & Uysdof_c = Uysdof[ii];

    h_cell_average   += hdof_c;
    n_cell_average   += ndof_c;
    Uwx_cell_average += Uxwdof_c;
    Uwy_cell_average += Uywdof_c;
    Usx_cell_average += Uxsdof_c;
    Usy_cell_average += Uysdof_c;
    
    
    fluxx_h_node[ii]   = h_flux_formula_x    (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_h_node[ii]   = h_flux_formula_y    (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_n_node[ii]   = n_flux_formula_x    (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_n_node[ii]   = n_flux_formula_y    (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Uwx_node[ii] = Uwx_flux_formula_x  (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Uwx_node[ii] = Uwx_flux_formula_y  (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Uwy_node[ii] = Uwy_flux_formula_x  (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Uwy_node[ii] = Uwy_flux_formula_y  (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Usx_node[ii] = Usx_flux_formula_x  (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Usx_node[ii] = Usx_flux_formula_y  (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Usy_node[ii] = Usy_flux_formula_x  (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Usy_node[ii] = Usy_flux_formula_y  (hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
  }
  h_cell_average     /= 4.;
  n_cell_average     /= 4.;
  Uwx_cell_average   /= 4.;
  Uwy_cell_average   /= 4.;
  Usx_cell_average   /= 4.;
  Usy_cell_average   /= 4.;
  

  
  const auto div_Fh_x = .5*((fluxx_h_node[1]-fluxx_h_node[0]) + (fluxx_h_node[3]-fluxx_h_node[2]));
  const auto div_Fh_y = .5*((fluxy_h_node[2]-fluxy_h_node[0]) + (fluxy_h_node[3]-fluxy_h_node[1]));
  const auto div_Fh_cell = Dy*div_Fh_x + Dx*div_Fh_y;

  const auto div_Fn_x = .5*((fluxx_n_node[1]-fluxx_n_node[0]) + (fluxx_n_node[3]-fluxx_n_node[2]));
  const auto div_Fn_y = .5*((fluxy_n_node[2]-fluxy_n_node[0]) + (fluxy_n_node[3]-fluxy_n_node[1]));
  const auto div_Fn_cell = Dy*div_Fn_x + Dx*div_Fn_y;
  
  const auto div_FUwx_x = .5*((fluxx_Uwx_node[1]-fluxx_Uwx_node[0]) + (fluxx_Uwx_node[3]-fluxx_Uwx_node[2]));
  const auto div_FUwx_y = .5*((fluxy_Uwx_node[2]-fluxy_Uwx_node[0]) + (fluxy_Uwx_node[3]-fluxy_Uwx_node[1]));
  const auto div_FUwx_cell = Dy*div_FUwx_x + Dx*div_FUwx_y;
  
  const auto div_FUwy_x = .5*((fluxx_Uwy_node[1]-fluxx_Uwy_node[0]) + (fluxx_Uwy_node[3]-fluxx_Uwy_node[2]));
  const auto div_FUwy_y = .5*((fluxy_Uwy_node[2]-fluxy_Uwy_node[0]) + (fluxy_Uwy_node[3]-fluxy_Uwy_node[1]));
  const auto div_FUwy_cell = Dy*div_FUwy_x + Dx*div_FUwy_y;

  const auto div_FUsx_x = .5*((fluxx_Usx_node[1]-fluxx_Usx_node[0]) + (fluxx_Usx_node[3]-fluxx_Usx_node[2]));
  const auto div_FUsx_y = .5*((fluxy_Usx_node[2]-fluxy_Usx_node[0]) + (fluxy_Usx_node[3]-fluxy_Usx_node[1]));
  const auto div_FUsx_cell = Dy*div_FUsx_x + Dx*div_FUsx_y;
  
  const auto div_FUsy_x = .5*((fluxx_Usy_node[1]-fluxx_Usy_node[0]) + (fluxx_Usy_node[3]-fluxx_Usy_node[2]));
  const auto div_FUsy_y = .5*((fluxy_Usy_node[2]-fluxy_Usy_node[0]) + (fluxy_Usy_node[3]-fluxy_Usy_node[1]));
  const auto div_FUsy_cell = Dy*div_FUsy_x + Dx*div_FUsy_y;


  const auto slope_x_c = ((Z_node[1] - Z_node[0]) + (Z_node[3] - Z_node[2]))/Dx/2.;
  const auto slope_y_c = ((Z_node[2] - Z_node[0]) + (Z_node[3] - Z_node[1]))/Dy/2.;

  const double divUw = ((Uxwdof[1] - Uxwdof[0]) + (Uxwdof[3] - Uxwdof[2]))/Dx/2. + ((Uywdof[2] - Uywdof[0]) + (Uywdof[3] - Uywdof[1]))/Dy/2.;
  const double divUs = ((Uxsdof[1] - Uxsdof[0]) + (Uxsdof[3] - Uxsdof[2]))/Dx/2. + ((Uysdof[2] - Uysdof[0]) + (Uysdof[3] - Uysdof[1]))/Dy/2.;

  const double tau = dt*.5;


  Z_onehalf[index_quadrant] = (Z_node[0]+Z_node[1]+Z_node[2]+Z_node[3])*.25;


  sol_onehalf[ordn     (index_quadrant)] = n_cell_average   - tau * (div_Fn_cell  /area + 1./(h_cell_average*(1.-n_cell_average))*divUw - 1./(h_cell_average*n_cell_average)*divUs);

  const auto v_h    = h_cell_average   - tau *  div_Fh_cell  /area;
  const auto v_Ux_w = Uwx_cell_average - tau * (div_FUwx_cell/area - src_slope_formula (    n_cell_average *h_cell_average, slope_x_c));
  const auto v_Uy_w = Uwy_cell_average - tau * (div_FUwy_cell/area - src_slope_formula (    n_cell_average *h_cell_average, slope_y_c));
  const auto v_Ux_s = Usx_cell_average - tau * (div_FUsx_cell/area - src_slope_formula ((1.-n_cell_average)*h_cell_average, slope_x_c));
  const auto v_Uy_s = Usy_cell_average - tau * (div_FUsy_cell/area - src_slope_formula ((1.-n_cell_average)*h_cell_average, slope_y_c));
  

  const auto & n_c = sol_onehalf[ordn     (index_quadrant)];

  const double density = (1.-n_c)*density_s + n_c*density_w;
  const double density_prime = (1.-n_c)*(density_s-density_w);


  const double tolerance = 1.e-4;
  const int Nmax = 1e3;
  int count = -1;
  double error = tolerance + 1; 
  while (count++<Nmax && error>tolerance)
  {
    const auto & h_c   = sol_onehalf[ordh    (index_quadrant)];
    const auto & Uxw_c = sol_onehalf[ordUxw  (index_quadrant)];
    const auto & Uyw_c = sol_onehalf[ordUyw  (index_quadrant)];
    const auto & Uxs_c = sol_onehalf[ordUxs  (index_quadrant)];
    const auto & Uys_c = sol_onehalf[ordUys  (index_quadrant)];

    const auto hw_c = n_c*h_c;
    const auto hs_c = (1.-n_c)*h_c; 


    const auto U_tot_x = Uxs_c + Uxw_c;
    const auto U_tot_y = Uys_c + Uyw_c;

    const auto abs_mass_flux = std::sqrt(U_tot_x*U_tot_x + U_tot_y*U_tot_y);

    const double Cw_d = hw_c>epsilon && h_c>epsilon ? (1.-n_c)/terminal_velocity/std::pow(n_c, m_coeff)*(density_s-density_w)*grav : 0.;
    const double Cs_d = hw_c>epsilon && h_c>epsilon ?     n_c /terminal_velocity/std::pow(n_c, m_coeff)*(density_s-density_w)*grav : 0.;


    const double vel_x_sign = std::abs(U_tot_x)>tolerance_sign ? U_tot_x/std::abs(U_tot_x) : U_tot_x/tolerance_sign;
    const double vel_y_sign = std::abs(U_tot_y)>tolerance_sign ? U_tot_y/std::abs(U_tot_y) : U_tot_y/tolerance_sign;


    const double Ctaux_der_h = (density_prime*grav*std::tan(bed_friction_angle_rad)*vel_x_sign - (h_c*h_c)>epsilon ? 2.*density*grav/turbulence_coeff/(h_c*h_c*h_c)*std::abs(U_tot_x)*U_tot_x : 0.);
    const double Ctauy_der_h = (density_prime*grav*std::tan(bed_friction_angle_rad)*vel_y_sign - (h_c*h_c)>epsilon ? 2.*density*grav/turbulence_coeff/(h_c*h_c*h_c)*std::abs(U_tot_y)*U_tot_y : 0.);

    // Newton Jacobian matrix, 
    // diagonalizzazione dei termini extra-diag delle matrici g e f
    const double a   = 1.;
    const double b_x = -tau*erosion_coefficient* (abs_mass_flux!=0 ? U_tot_x/abs_mass_flux : 0.);
    const double b_y = -tau*erosion_coefficient* (abs_mass_flux!=0 ? U_tot_y/abs_mass_flux : 0.);
    const double c_x = b_x;
    const double c_y = b_y;
    const double d   = 1. + tau* Cw_d       /density_w;
    const double e   =    - tau* Cs_d       /density_w;
    const double f_x =      tau* Ctaux_der_h/density_s;
    const double f_y =      tau* Ctauy_der_h/density_s;
    const double g_x =      tau* (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign) + (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_x) : 0.) /density_s - tau*Cw_d/density_s;
    const double g_y =      tau* (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign) + (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_y) : 0.) /density_s - tau*Cw_d/density_s;
    const double h_x = 1. - tau* (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign) + (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_x) : 0.) /density_s + tau*Cs_d/density_s;
    const double h_y = 1. - tau* (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign) + (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_y) : 0.) /density_s + tau*Cs_d/density_s;

    // rhs terms of the Newton matrix,
    const auto f_h    = v_h    + tau*h_src_formula  (h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - h_c;
    const auto f_Ux_w = v_Ux_w + tau*Uxw_src_formula(h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - Uxw_c;
    const auto f_Uy_w = v_Uy_w + tau*Uyw_src_formula(h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - Uyw_c;
    const auto f_Ux_s = v_Ux_s + tau*Uxs_src_formula(h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - Uxs_c;
    const auto f_Uy_s = v_Uy_s + tau*Uys_src_formula(h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - Uys_c;

    // solve the lin. system,
    const double big_A = -f_x*c_x/a + f_x*b_x/d/a*e - g_x/d*e + h_x;
    const double big_B = -f_x*c_y/a + f_x*b_y/d/a*e;
    const double big_C = -f_y*c_x/a + f_y*b_x/d/a*e;
    const double big_D = -f_y*c_y/a + f_y*b_y/d/a*e - g_y/d*e + h_y;

    const auto rhs_x = f_Ux_s - f_x*f_h/a + (f_x*b_x*f_Ux_w+f_x*b_y*f_Uy_w)/d/a - g_x/d*f_Ux_w;
    const auto rhs_y = f_Uy_s - f_y*f_h/a + (f_y*b_x*f_Ux_w+f_y*b_y*f_Uy_w)/d/a - g_y/d*f_Uy_w;

    const double big_det = big_A*big_D-big_B*big_C;

    const auto delta_Ux_s = (big_D*rhs_x-big_B*rhs_y)/big_det;
    const auto delta_Uy_s = (big_A*rhs_y-big_C*rhs_x)/big_det; 

    const auto delta_Ux_w = (f_Ux_w-e*delta_Ux_s)/d;
    const auto delta_Uy_w = (f_Uy_w-e*delta_Uy_s)/d;

    const auto delta_h    = (f_h - (b_x*delta_Ux_w+b_y*delta_Uy_w) - (c_x*delta_Ux_s+c_y*delta_Uy_s))/a;
    

    error = std::sqrt(delta_h*delta_h + delta_Ux_s*delta_Ux_s + delta_Uy_s*delta_Uy_s + delta_Ux_w*delta_Ux_w + delta_Uy_w*delta_Uy_w);

    sol_onehalf[ordh    (index_quadrant)] += delta_h;
    sol_onehalf[ordUxw  (index_quadrant)] += delta_Ux_w;
    sol_onehalf[ordUyw  (index_quadrant)] += delta_Uy_w;
    sol_onehalf[ordUxs  (index_quadrant)] += delta_Ux_s;
    sol_onehalf[ordUys  (index_quadrant)] += delta_Uy_s;

  }

  // apply the corrector step now,
  const double r_coeff = density_w/density_s;
  const double g_prime = (1.-r_coeff)*grav;

  const auto & h_c   = sol_onehalf[ordh    (index_quadrant)];
  const auto & Uxw_c = sol_onehalf[ordUxw  (index_quadrant)];
  const auto & Uyw_c = sol_onehalf[ordUyw  (index_quadrant)];
  const auto & Uxs_c = sol_onehalf[ordUxs  (index_quadrant)];
  const auto & Uys_c = sol_onehalf[ordUys  (index_quadrant)];

  const double hw_c = n_c*h_c;
  const double hs_c = (1.-n_c)*h_c; 

  const double vel_w_x = hw_c>epsilon ? Uxw_c/hw_c : 0.;
  const double vel_s_x = hs_c>epsilon ? Uxs_c/hs_c : 0.;
  const double vel_w_y = hw_c>epsilon ? Uyw_c/hw_c : 0.;
  const double vel_s_y = hs_c>epsilon ? Uys_c/hs_c : 0.;

  const double cx_sgn = hw_c*hs_c/(tau*(hs_c+r_coeff*hw_c))* std::max(std::abs(vel_w_x-vel_s_x)/std::sqrt(g_prime*(hw_c+hs_c)) - 1., 0.);
  const double cy_sgn = hw_c*hs_c/(tau*(hs_c+r_coeff*hw_c))* std::max(std::abs(vel_w_y-vel_s_y)/std::sqrt(g_prime*(hw_c+hs_c)) - 1., 0.);

  const double big_Ax = 1. + hw_c > epsilon ? tau*cx_sgn/hw_c : 0.;
  const double big_Bx = hs_c > epsilon ? -tau*cx_sgn/hs_c : 0.;
  const double big_Cx = hw_c > epsilon ? -tau*r_coeff*cx_sgn/hw_c : 0.;
  const double big_Dx = 1. + hs_c > epsilon ? tau*r_coeff*cx_sgn/hs_c : 0.;

  const double rhs_1x = Uxw_c;
  const double rhs_2x = Uxs_c;

  const double big_detx = big_Ax*big_Dx-big_Bx*big_Cx;

  const double big_Ay = 1. + hw_c > epsilon ? tau*cy_sgn/hw_c : 0.;
  const double big_By = hs_c > epsilon ? -tau*cy_sgn/hs_c : 0.;
  const double big_Cy = hw_c > epsilon ? -tau*r_coeff*cy_sgn/hw_c : 0.;
  const double big_Dy = 1. + hs_c > epsilon ? tau*r_coeff*cy_sgn/hs_c : 0.;

  const double rhs_1y = Uyw_c;
  const double rhs_2y = Uys_c;

  const double big_dety = big_Ay*big_Dy-big_By*big_Cy;


  sol_onehalf[ordUxw  (index_quadrant)] = (rhs_1x*big_Dx-rhs_2x*big_Bx)/big_detx;
  sol_onehalf[ordUxs  (index_quadrant)] = (rhs_2x*big_Ax-rhs_1x*big_Cx)/big_detx;

  sol_onehalf[ordUyw  (index_quadrant)] = (rhs_1y*big_Dy-rhs_2y*big_By)/big_dety;
  sol_onehalf[ordUys  (index_quadrant)] = (rhs_2y*big_Ay-rhs_1y*big_Cy)/big_dety;
  
  

}



void
TG2_scheme::solve_non_lin(const int& kk)
{
  const auto & hold_c   = sold.get_owned_data ()[kk  ];
  const auto & nold_c   = sold.get_owned_data ()[kk+1];
  const auto & Uxwold_c = sold.get_owned_data ()[kk+2];
  const auto & Uywold_c = sold.get_owned_data ()[kk+3];
  const auto & Uxsold_c = sold.get_owned_data ()[kk+4];
  const auto & Uysold_c = sold.get_owned_data ()[kk+5];


  sol.get_owned_data ()[kk+1] += dt*incr.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1];


  // solve non-linearities like the first step of the TG2 method to get the complete low order solution,
  // terminated this part add the nodal correction to get the hyperbolicity,
  const auto v_h    = sol.get_owned_data ()[kk  ] + dt*incr.get_owned_data ()[kk  ]/mass.get_owned_data ()[kk  ] + dt*.5*h_src_formula  (hold_c, nold_c, Uxwold_c, Uywold_c, Uxsold_c, Uysold_c);
  const auto v_Ux_w = sol.get_owned_data ()[kk+2] + dt*incr.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] + dt*.5*Uxw_src_formula(hold_c, nold_c, Uxwold_c, Uywold_c, Uxsold_c, Uysold_c);
  const auto v_Uy_w = sol.get_owned_data ()[kk+3] + dt*incr.get_owned_data ()[kk+3]/mass.get_owned_data ()[kk+3] + dt*.5*Uyw_src_formula(hold_c, nold_c, Uxwold_c, Uywold_c, Uxsold_c, Uysold_c);
  const auto v_Ux_s = sol.get_owned_data ()[kk+4] + dt*incr.get_owned_data ()[kk+4]/mass.get_owned_data ()[kk+4] + dt*.5*Uxs_src_formula(hold_c, nold_c, Uxwold_c, Uywold_c, Uxsold_c, Uysold_c);
  const auto v_Uy_s = sol.get_owned_data ()[kk+5] + dt*incr.get_owned_data ()[kk+5]/mass.get_owned_data ()[kk+5] + dt*.5*Uys_src_formula(hold_c, nold_c, Uxwold_c, Uywold_c, Uxsold_c, Uysold_c);


  const auto & n_c = sol.get_owned_data ()[kk+1];

  const double density = (1.-n_c)*density_s + n_c*density_w;
  const double density_prime = (1.-n_c)*(density_s-density_w);

  double tau = dt*.5;


  const double tolerance = 1.e-4;
  const int Nmax = 1e3;
  int count = -1;
  double error = tolerance + 1;
  while (count++<Nmax && error>tolerance)
  {
    const auto & h_c   = sol.get_owned_data ()[kk  ];
    const auto & Uxw_c = sol.get_owned_data ()[kk+2];
    const auto & Uyw_c = sol.get_owned_data ()[kk+3];
    const auto & Uxs_c = sol.get_owned_data ()[kk+4];
    const auto & Uys_c = sol.get_owned_data ()[kk+5];

    const auto hw_c = n_c*h_c;
    const auto hs_c = (1.-n_c)*h_c; 


    const auto U_tot_x = Uxs_c + Uxw_c;
    const auto U_tot_y = Uys_c + Uyw_c;

    const auto abs_mass_flux = std::sqrt(U_tot_x*U_tot_x + U_tot_y*U_tot_y);

    const double Cw_d = hw_c>epsilon && h_c>epsilon ? (1.-n_c)/terminal_velocity/std::pow(n_c, m_coeff)*(density_s-density_w)*grav : 0.;
    const double Cs_d = hw_c>epsilon && h_c>epsilon ?     n_c /terminal_velocity/std::pow(n_c, m_coeff)*(density_s-density_w)*grav : 0.;


    const double vel_x_sign = std::abs(U_tot_x)>tolerance_sign ? U_tot_x/std::abs(U_tot_x) : U_tot_x/tolerance_sign;
    const double vel_y_sign = std::abs(U_tot_y)>tolerance_sign ? U_tot_y/std::abs(U_tot_y) : U_tot_y/tolerance_sign;


    const double Ctaux_der_h = (density_prime*grav*std::tan(bed_friction_angle_rad)*vel_x_sign - (h_c*h_c)>epsilon ? 2.*density*grav/turbulence_coeff/(h_c*h_c*h_c)*std::abs(U_tot_x)*U_tot_x : 0.);
    const double Ctauy_der_h = (density_prime*grav*std::tan(bed_friction_angle_rad)*vel_y_sign - (h_c*h_c)>epsilon ? 2.*density*grav/turbulence_coeff/(h_c*h_c*h_c)*std::abs(U_tot_y)*U_tot_y : 0.);

    // Newton Jacobian matrix, 
    // diagonalizzazione dei termini extra-diag delle matrici g e f
    const double a   = 1.;
    const double b_x = -tau*erosion_coefficient* (abs_mass_flux!=0 ? U_tot_x/abs_mass_flux : 0.);
    const double b_y = -tau*erosion_coefficient* (abs_mass_flux!=0 ? U_tot_y/abs_mass_flux : 0.);
    const double c_x = b_x;
    const double c_y = b_y;
    const double d   = 1. + tau* Cw_d       /density_w;
    const double e   =    - tau* Cs_d       /density_w;
    const double f_x =      tau* Ctaux_der_h/density_s;
    const double f_y =      tau* Ctauy_der_h/density_s;
    const double g_x =      tau* (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign) + (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_x) : 0.) /density_s - tau*Cw_d/density_s;
    const double g_y =      tau* (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign) + (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_y) : 0.) /density_s - tau*Cw_d/density_s;
    const double h_x = 1. - tau* (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign) + (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_x) : 0.) /density_s + tau*Cs_d/density_s;
    const double h_y = 1. - tau* (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*(std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign) + (h_c*h_c)>epsilon ? density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_y) : 0.) /density_s + tau*Cs_d/density_s;

    // rhs terms of the Newton matrix,
    const auto f_h    = v_h    + tau*h_src_formula  (h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - h_c;
    const auto f_Ux_w = v_Ux_w + tau*Uxw_src_formula(h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - Uxw_c;
    const auto f_Uy_w = v_Uy_w + tau*Uyw_src_formula(h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - Uyw_c;
    const auto f_Ux_s = v_Ux_s + tau*Uxs_src_formula(h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - Uxs_c;
    const auto f_Uy_s = v_Uy_s + tau*Uys_src_formula(h_c, n_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - Uys_c;

    // solve the lin. system,
    const double big_A = -f_x*c_x/a + f_x*b_x/d/a*e - g_x/d*e + h_x;
    const double big_B = -f_x*c_y/a + f_x*b_y/d/a*e;
    const double big_C = -f_y*c_x/a + f_y*b_x/d/a*e;
    const double big_D = -f_y*c_y/a + f_y*b_y/d/a*e - g_y/d*e + h_y;

    const auto rhs_x = f_Ux_s - f_x*f_h/a + (f_x*b_x*f_Ux_w+f_x*b_y*f_Uy_w)/d/a - g_x/d*f_Ux_w;
    const auto rhs_y = f_Uy_s - f_y*f_h/a + (f_y*b_x*f_Ux_w+f_y*b_y*f_Uy_w)/d/a - g_y/d*f_Uy_w;

    const double big_det = big_A*big_D-big_B*big_C;

    const auto delta_Ux_s = (big_D*rhs_x-big_B*rhs_y)/big_det;
    const auto delta_Uy_s = (big_A*rhs_y-big_C*rhs_x)/big_det;

    const auto delta_Ux_w = (f_Ux_w-e*delta_Ux_s)/d;
    const auto delta_Uy_w = (f_Uy_w-e*delta_Uy_s)/d;

    const auto delta_h    = (f_h - (b_x*delta_Ux_w+b_y*delta_Uy_w) - (c_x*delta_Ux_s+c_y*delta_Uy_s))/a;
    
    

    error = std::sqrt(delta_h*delta_h + delta_Ux_s*delta_Ux_s + delta_Uy_s*delta_Uy_s + delta_Ux_w*delta_Ux_w + delta_Uy_w*delta_Uy_w);

    sol.get_owned_data ()[kk  ] += delta_h;
    sol.get_owned_data ()[kk+2] += delta_Ux_w;
    sol.get_owned_data ()[kk+3] += delta_Uy_w;
    sol.get_owned_data ()[kk+4] += delta_Ux_s;
    sol.get_owned_data ()[kk+5] += delta_Uy_s;
  }


  // apply the corrector step now, in this way we obtain a complete low order solution, hopefully we are killing the oscillations of the first order method,
  tau = dt;

  const double r_coeff = density_w/density_s;
  const double g_prime = (1.-r_coeff)*grav;

  const auto & h_c   = sol.get_owned_data ()[kk  ];
  const auto & Uxw_c = sol.get_owned_data ()[kk+2];
  const auto & Uyw_c = sol.get_owned_data ()[kk+3];
  const auto & Uxs_c = sol.get_owned_data ()[kk+4];
  const auto & Uys_c = sol.get_owned_data ()[kk+5];

  const double hw_c = n_c*h_c;
  const double hs_c = (1.-n_c)*h_c; 

  const double vel_w_x = hw_c>epsilon ? Uxw_c/hw_c : 0.;
  const double vel_s_x = hs_c>epsilon ? Uxs_c/hs_c : 0.;
  const double vel_w_y = hw_c>epsilon ? Uyw_c/hw_c : 0.;
  const double vel_s_y = hs_c>epsilon ? Uys_c/hs_c : 0.;

  const double cx_sgn = hw_c*hs_c/(tau*(hs_c+r_coeff*hw_c))* std::max(std::abs(vel_w_x-vel_s_x)/std::sqrt(g_prime*(hw_c+hs_c)) - 1., 0.);
  const double cy_sgn = hw_c*hs_c/(tau*(hs_c+r_coeff*hw_c))* std::max(std::abs(vel_w_y-vel_s_y)/std::sqrt(g_prime*(hw_c+hs_c)) - 1., 0.);

  const double big_Ax = 1. + hw_c > epsilon ? tau*cx_sgn/hw_c : 0.;
  const double big_Bx = hs_c > epsilon ? -tau*cx_sgn/hs_c : 0.;
  const double big_Cx = hw_c > epsilon ? -tau*r_coeff*cx_sgn/hw_c : 0.;
  const double big_Dx = 1. + hs_c > epsilon ? tau*r_coeff*cx_sgn/hs_c : 0.;

  const double rhs_1x = Uxw_c;
  const double rhs_2x = Uxs_c;

  const double big_detx = big_Ax*big_Dx-big_Bx*big_Cx;

  const double big_Ay = 1. + hw_c > epsilon ? tau*cy_sgn/hw_c : 0.;
  const double big_By = hs_c > epsilon ? -tau*cy_sgn/hs_c : 0.;
  const double big_Cy = hw_c > epsilon ? -tau*r_coeff*cy_sgn/hw_c : 0.;
  const double big_Dy = 1. + hs_c > epsilon ? tau*r_coeff*cy_sgn/hs_c : 0.;

  const double rhs_1y = Uyw_c;
  const double rhs_2y = Uys_c;

  const double big_dety = big_Ay*big_Dy-big_By*big_Cy;


  sol.get_owned_data ()[kk+2] = (rhs_1x*big_Dx-rhs_2x*big_Bx)/big_detx;
  sol.get_owned_data ()[kk+3] = (rhs_1y*big_Dy-rhs_2y*big_By)/big_dety;
  sol.get_owned_data ()[kk+4] = (rhs_2x*big_Ax-rhs_1x*big_Cx)/big_detx;
  sol.get_owned_data ()[kk+5] = (rhs_2y*big_Ay-rhs_1y*big_Cy)/big_dety;



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

  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hdof[ii]   = sol [ordh   (quadrant->gt (ii) )];
      ndof[ii]   = sol [ordn   (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 

      Z_node[ii]  = Z [quadrant->gt (ii)];

      isdof_or_hanging[ii] = 1.;
    }
    else
    {
      hdof[ii]   = .5 * (sol [ordh   (quadrant->gparent (0, ii) )] +
                         sol [ordh   (quadrant->gparent (1, ii) )]);
      ndof[ii]   = .5 * (sol [ordn   (quadrant->gparent (0, ii) )] +
                         sol [ordn   (quadrant->gparent (1, ii) )]);
      Uxwdof[ii] = .5 * (sol [ordUxw (quadrant->gparent (0, ii) )] +
                         sol [ordUxw (quadrant->gparent (1, ii) )]);
      Uywdof[ii] = .5 * (sol [ordUyw (quadrant->gparent (0, ii) )] +
                         sol [ordUyw (quadrant->gparent (1, ii) )]);
      Uxsdof[ii] = .5 * (sol [ordUxs (quadrant->gparent (0, ii) )] +
                         sol [ordUxs (quadrant->gparent (1, ii) )]);
      Uysdof[ii] = .5 * (sol [ordUys (quadrant->gparent (0, ii) )] +
                         sol [ordUys (quadrant->gparent (1, ii) )]);

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
    const auto lambdas = max_eigen (hdof[ii], ndof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii]);
    vel_rusanov_cell_x += lambdas[0];
    vel_rusanov_cell_y += lambdas[1];
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.; 

  const std::array<double,4> eta_vec = {hdof[0]+Z_node[0], hdof[1]+Z_node[1], hdof[2]+Z_node[2], hdof[3]+Z_node[3]};
 
  grad_cell_eta   = {.5 * ( (eta_vec [3] - eta_vec [2]) + (eta_vec [1] - eta_vec [0]) ), .5 * ( (eta_vec [2] - eta_vec [0]) + (eta_vec [3] - eta_vec [1]) )};
  grad_cell_h     = {.5 * ( (hdof    [3] - hdof    [2]) + (hdof    [1] - hdof    [0]) ), .5 * ( (hdof    [2] - hdof    [0]) + (hdof    [3] - hdof    [1]) )};
  grad_cell_n     = {.5 * ( (ndof    [3] - ndof    [2]) + (ndof    [1] - ndof    [0]) ), .5 * ( (ndof    [2] - ndof    [0]) + (ndof    [3] - ndof    [1]) )};
  grad_cell_Uxw   = {.5 * ( (Uxwdof  [3] - Uxwdof  [2]) + (Uxwdof  [1] - Uxwdof  [0]) ), .5 * ( (Uxwdof  [2] - Uxwdof  [0]) + (Uxwdof  [3] - Uxwdof  [1]) )};
  grad_cell_Uyw   = {.5 * ( (Uywdof  [3] - Uywdof  [2]) + (Uywdof  [1] - Uywdof  [0]) ), .5 * ( (Uywdof  [2] - Uywdof  [0]) + (Uywdof  [3] - Uywdof  [1]) )};
  grad_cell_Uxs   = {.5 * ( (Uxsdof  [3] - Uxsdof  [2]) + (Uxsdof  [1] - Uxsdof  [0]) ), .5 * ( (Uxsdof  [2] - Uxsdof  [0]) + (Uxsdof  [3] - Uxsdof  [1]) )};
  grad_cell_Uys   = {.5 * ( (Uysdof  [3] - Uysdof  [2]) + (Uysdof  [1] - Uysdof  [0]) ), .5 * ( (Uysdof  [2] - Uysdof  [0]) + (Uysdof  [3] - Uysdof  [1]) )};
 
  //std::cout << grad_cell_eta[0] << " " << grad_cell_eta[1] << " " << grad_cell_eta[2] << " " << grad_cell_eta[3] << std::endl;

  const double & h_cell     = sol_onehalf[ordh     (index_quadrant)];
  const double & n_cell     = sol_onehalf[ordn     (index_quadrant)];
  const double & Uxw_cell   = sol_onehalf[ordUxw   (index_quadrant)];
  const double & Uyw_cell   = sol_onehalf[ordUyw   (index_quadrant)];
  const double & Uxs_cell   = sol_onehalf[ordUxs   (index_quadrant)];
  const double & Uys_cell   = sol_onehalf[ordUys   (index_quadrant)];

  const double & Z_cell     = Z_onehalf[index_quadrant];

  const auto hs_cell = (1.-n_cell)*h_cell;
  const auto hw_cell = n_cell*h_cell;

  std::array<double,4> contr_x_n = {0., 0., 0., 0.}, contr_y_n = {0., 0., 0., 0.}, contr_x_w = {0., 0., 0., 0.}, contr_y_w = {0., 0., 0., 0.}, contr_x_s = {0., 0., 0., 0.}, contr_y_s = {0., 0., 0., 0.};

  const auto dens_ratio = density_w/density_s;


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
          auto Z_cell_nei = Z_onehalf[index_quadrant_nei];

          auto h_cell_nei     = sol_onehalf[ordh     (index_quadrant)];
          auto n_cell_nei     = sol_onehalf[ordn     (index_quadrant)];
          auto Uxw_cell_nei   = sol_onehalf[ordUxw   (index_quadrant)];
          auto Uyw_cell_nei   = sol_onehalf[ordUyw   (index_quadrant)];
          auto Uxs_cell_nei   = sol_onehalf[ordUxs   (index_quadrant)];
          auto Uys_cell_nei   = sol_onehalf[ordUys   (index_quadrant)];


          // .5 salta fuori dall'integrazione per trapezi tra 0 e 1 in coordinata \xi (è il valore in LHS da metter qui sotto!)
          contr_x_n[i_1] += .5*signum(outward_normal_edge[0])*((hs_cell>epsilon ? 1./hs_cell : 0.)*(Uxw_cell_nei - Uxw_cell) + (hw_cell>epsilon ? -1./hw_cell : 0.)*(Uxs_cell_nei - Uxs_cell))*isdof_or_hanging[i_1];
          contr_x_n[i_2] += .5*signum(outward_normal_edge[0])*((hs_cell>epsilon ? 1./hs_cell : 0.)*(Uxw_cell_nei - Uxw_cell) + (hw_cell>epsilon ? -1./hw_cell : 0.)*(Uxs_cell_nei - Uxs_cell))*isdof_or_hanging[i_2];

          contr_y_n[i_1] += .5*signum(outward_normal_edge[1])*((hs_cell>epsilon ? 1./hs_cell : 0.)*(Uxw_cell_nei - Uxw_cell) + (hw_cell>epsilon ? -1./hw_cell : 0.)*(Uxs_cell_nei - Uxs_cell))*isdof_or_hanging[i_1];
          contr_y_n[i_2] += .5*signum(outward_normal_edge[1])*((hs_cell>epsilon ? 1./hs_cell : 0.)*(Uxw_cell_nei - Uxw_cell) + (hw_cell>epsilon ? -1./hw_cell : 0.)*(Uxs_cell_nei - Uxs_cell))*isdof_or_hanging[i_2];

          contr_x_w[i_1] += .5*signum(outward_normal_edge[0])*(grav*n_cell*h_cell*(Z_cell_nei - Z_cell) -.5*grav*h_cell*h_cell*(n_cell_nei - n_cell) )*isdof_or_hanging[i_1];
          contr_x_w[i_2] += .5*signum(outward_normal_edge[0])*(grav*n_cell*h_cell*(Z_cell_nei - Z_cell) -.5*grav*h_cell*h_cell*(n_cell_nei - n_cell) )*isdof_or_hanging[i_2];

          contr_y_w[i_1] += .5*signum(outward_normal_edge[1])*(grav*n_cell*h_cell*(Z_cell_nei - Z_cell) -.5*grav*h_cell*h_cell*(n_cell_nei - n_cell) )*isdof_or_hanging[i_1];
          contr_y_w[i_2] += .5*signum(outward_normal_edge[1])*(grav*n_cell*h_cell*(Z_cell_nei - Z_cell) -.5*grav*h_cell*h_cell*(n_cell_nei - n_cell) )*isdof_or_hanging[i_2];

          contr_x_s[i_1] += .5*signum(outward_normal_edge[0])*(grav*(1.-n_cell)*h_cell*(Z_cell_nei - Z_cell) +.5*grav*h_cell*h_cell*(n_cell_nei - n_cell)*dens_ratio )*isdof_or_hanging[i_1];
          contr_x_s[i_2] += .5*signum(outward_normal_edge[0])*(grav*(1.-n_cell)*h_cell*(Z_cell_nei - Z_cell) +.5*grav*h_cell*h_cell*(n_cell_nei - n_cell)*dens_ratio )*isdof_or_hanging[i_2];

          contr_y_s[i_1] += .5*signum(outward_normal_edge[1])*(grav*(1.-n_cell)*h_cell*(Z_cell_nei - Z_cell) +.5*grav*h_cell*h_cell*(n_cell_nei - n_cell)*dens_ratio )*isdof_or_hanging[i_1];
          contr_y_s[i_2] += .5*signum(outward_normal_edge[1])*(grav*(1.-n_cell)*h_cell*(Z_cell_nei - Z_cell) +.5*grav*h_cell*h_cell*(n_cell_nei - n_cell)*dens_ratio )*isdof_or_hanging[i_2];

          //break; // this just goes outside the jEdge cycle 
        }
      }

    }

    if (is_boundary_edge) // set boundary conditions
    { 
      
      auto h_cell_nei   = h_cell;
      auto n_cell_nei   = n_cell;
      auto Uxw_cell_nei = Uxw_cell;
      auto Uyw_cell_nei = Uyw_cell;
      auto Uxs_cell_nei = Uxs_cell;
      auto Uys_cell_nei = Uys_cell;

      Uxw_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxw_cell + outward_normal_edge[1]*Uyw_cell)*outward_normal_edge[0];
      Uyw_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxw_cell + outward_normal_edge[1]*Uyw_cell)*outward_normal_edge[1];

      Uxs_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxs_cell + outward_normal_edge[1]*Uys_cell)*outward_normal_edge[0];
      Uys_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxs_cell + outward_normal_edge[1]*Uys_cell)*outward_normal_edge[1];

      const auto speed     = max_eigen(h_cell,     n_cell,     Uxw_cell,     Uyw_cell,     Uxs_cell,     Uys_cell    );
      const auto speed_nei = max_eigen(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei);

      const auto smax = std::max(speed[0]*outward_normal_edge[0]+speed[1]*outward_normal_edge[1], speed_nei[0]*outward_normal_edge[0]+speed_nei[1]*outward_normal_edge[1]); 

      const auto flux_int_h   = .5*((h_flux_formula_x  (h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+h_flux_formula_x  (h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (h_flux_formula_y  (h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+h_flux_formula_y  (h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(h_cell_nei  -h_cell  );
      const auto flux_int_n   = .5*((n_flux_formula_x  (h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+n_flux_formula_x  (h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (n_flux_formula_y  (h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+n_flux_formula_y  (h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(n_cell_nei  -n_cell  );
      const auto flux_int_Uxw = .5*((Uwx_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uwx_flux_formula_x(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (Uwx_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uwx_flux_formula_y(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uxw_cell_nei-Uxw_cell);
      const auto flux_int_Uyw = .5*((Uwy_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uwy_flux_formula_x(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (Uwy_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uwy_flux_formula_y(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uyw_cell_nei-Uyw_cell);
      const auto flux_int_Uxs = .5*((Usx_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Usx_flux_formula_x(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (Usx_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Usx_flux_formula_y(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uxs_cell_nei-Uxs_cell);
      const auto flux_int_Uys = .5*((Usy_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Usy_flux_formula_x(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (Usy_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Usy_flux_formula_y(h_cell_nei, n_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uys_cell_nei-Uys_cell);


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
        incr[ordh    (quadrant->gt (i_1))] += -edge_length*flux_int_h  *.5;
        incr[ordn    (quadrant->gt (i_1))] += -edge_length*flux_int_n  *.5;
        incr[ordUxw  (quadrant->gt (i_1))] += -edge_length*flux_int_Uxw*.5;
        incr[ordUyw  (quadrant->gt (i_1))] += -edge_length*flux_int_Uyw*.5;
        incr[ordUxs  (quadrant->gt (i_1))] += -edge_length*flux_int_Uxs*.5;
        incr[ordUys  (quadrant->gt (i_1))] += -edge_length*flux_int_Uys*.5;
      }
      else
      {
        incr [ordh   (quadrant->gparent(0,i_1))] += -edge_length*flux_int_h  *.5*.5;
        incr [ordh   (quadrant->gparent(1,i_1))] += -edge_length*flux_int_h  *.5*.5;

        incr [ordn   (quadrant->gparent(0,i_1))] += -edge_length*flux_int_n  *.5*.5;
        incr [ordn   (quadrant->gparent(1,i_1))] += -edge_length*flux_int_n  *.5*.5;
      
        incr [ordUxw (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uxw*.5*.5;
        incr [ordUxw (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uxw*.5*.5;
      
        incr [ordUyw (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uyw*.5*.5;
        incr [ordUyw (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uyw*.5*.5;

        incr [ordUxs (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uxs*.5*.5;
        incr [ordUxs (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uxs*.5*.5;
      
        incr [ordUys (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uys*.5*.5;
        incr [ordUys (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uys*.5*.5;
      }

      if (! quadrant->is_hanging (i_2))
      {
        incr[ordh    (quadrant->gt (i_2))] += -edge_length*flux_int_h  *.5;
        incr[ordn    (quadrant->gt (i_2))] += -edge_length*flux_int_n  *.5;
        incr[ordUxw  (quadrant->gt (i_2))] += -edge_length*flux_int_Uxw*.5;
        incr[ordUyw  (quadrant->gt (i_2))] += -edge_length*flux_int_Uyw*.5;
        incr[ordUxs  (quadrant->gt (i_2))] += -edge_length*flux_int_Uxs*.5;
        incr[ordUys  (quadrant->gt (i_2))] += -edge_length*flux_int_Uys*.5;
      }
      else
      {
        // il secondo .5 è per hanging nodes
        incr [ordh   (quadrant->gparent(0,i_2))] += -edge_length*flux_int_h  *.5*.5;
        incr [ordh   (quadrant->gparent(1,i_2))] += -edge_length*flux_int_h  *.5*.5;

        incr [ordn   (quadrant->gparent(0,i_2))] += -edge_length*flux_int_n  *.5*.5;
        incr [ordn   (quadrant->gparent(1,i_2))] += -edge_length*flux_int_n  *.5*.5;
      
        incr [ordUxw (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uxw*.5*.5;
        incr [ordUxw (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uxw*.5*.5;
      
        incr [ordUyw (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uyw*.5*.5;
        incr [ordUyw (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uyw*.5*.5;

        incr [ordUxs (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uxs*.5*.5;
        incr [ordUxs (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uxs*.5*.5;
      
        incr [ordUys (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uys*.5*.5;
        incr [ordUys (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uys*.5*.5;
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
  

  const auto diff_term_h_x   = h_cell>epsilon ? grad_cell_eta[0]*vel_rusanov_cell_y*.5 : grad_cell_h[0]*vel_rusanov_cell_y*.5;
  const auto diff_term_h_y   = h_cell>epsilon ? grad_cell_eta[1]*vel_rusanov_cell_x*.5 : grad_cell_h[1]*vel_rusanov_cell_x*.5;

  const auto diff_term_n_x   = grad_cell_n [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_n_y   = grad_cell_n [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uxw_x = grad_cell_Uxw [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uxw_y = grad_cell_Uxw [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uyw_x = grad_cell_Uyw [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uyw_y = grad_cell_Uyw [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uxs_x = grad_cell_Uxs [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uxs_y = grad_cell_Uxs [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uys_x = grad_cell_Uys [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uys_y = grad_cell_Uys [1]*vel_rusanov_cell_x*.5;




  const auto F_star_h_x   = h_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_h_x;
  const auto F_star_h_y   = h_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_h_y;

  const auto F_star_n_x   = n_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_n_x;
  const auto F_star_n_y   = n_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_n_y;

  const auto F_star_Uxw_x = Uwx_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uxw_x;
  const auto F_star_Uxw_y = Uwx_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uxw_y;

  const auto F_star_Uyw_x = Uwy_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uyw_x;
  const auto F_star_Uyw_y = Uwy_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uyw_y;

  const auto F_star_Uxs_x = Usx_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uxs_x;
  const auto F_star_Uxs_y = Usx_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uxs_y;

  const auto F_star_Uys_x = Usy_flux_formula_x(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uys_x;
  const auto F_star_Uys_y = Usy_flux_formula_y(h_cell, n_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uys_y;


  for (int ii = 0; ii < 4; ++ii){

    const auto h_     = der_coeffs_x[ii]*F_star_h_x  +der_coeffs_y[ii]*F_star_h_y;
    const auto n_     = der_coeffs_x[ii]*F_star_n_x  +der_coeffs_y[ii]*F_star_n_y   - .5*Dy*contr_x_n[ii] - .5*Dx*contr_y_n[ii];
    const auto Uxw_   = der_coeffs_x[ii]*F_star_Uxw_x+der_coeffs_y[ii]*F_star_Uxw_y - .5*Dy*contr_x_w[ii];
    const auto Uyw_   = der_coeffs_x[ii]*F_star_Uyw_x+der_coeffs_y[ii]*F_star_Uyw_y - .5*Dx*contr_y_w[ii];
    const auto Uxs_   = der_coeffs_x[ii]*F_star_Uxs_x+der_coeffs_y[ii]*F_star_Uxs_y - .5*Dy*contr_x_s[ii];
    const auto Uys_   = der_coeffs_x[ii]*F_star_Uys_x+der_coeffs_y[ii]*F_star_Uys_y - .5*Dx*contr_y_s[ii];

    const auto h_al   = der_coeffs_x[ii]*diff_term_h_x   + der_coeffs_y[ii]*diff_term_h_y; 
    const auto n_al   = der_coeffs_x[ii]*diff_term_n_x   + der_coeffs_y[ii]*diff_term_n_y; 
    const auto Uxw_al = der_coeffs_x[ii]*diff_term_Uxw_x + der_coeffs_y[ii]*diff_term_Uxw_y;
    const auto Uyw_al = der_coeffs_x[ii]*diff_term_Uyw_x + der_coeffs_y[ii]*diff_term_Uyw_y;
    const auto Uxs_al = der_coeffs_x[ii]*diff_term_Uxs_x + der_coeffs_y[ii]*diff_term_Uxs_y;
    const auto Uys_al = der_coeffs_x[ii]*diff_term_Uys_x + der_coeffs_y[ii]*diff_term_Uys_y;

    incr_anti_diff[ordh  (index_quadrant)][ii] = h_al;
    incr_anti_diff[ordn  (index_quadrant)][ii] = n_al;
    incr_anti_diff[ordUxw(index_quadrant)][ii] = Uxw_al;
    incr_anti_diff[ordUyw(index_quadrant)][ii] = Uyw_al;
    incr_anti_diff[ordUxs(index_quadrant)][ii] = Uxs_al;
    incr_anti_diff[ordUys(index_quadrant)][ii] = Uys_al;


    if (! quadrant->is_hanging (ii)){

      incr [ordh   (quadrant->gt (ii))] += h_;
      incr [ordn   (quadrant->gt (ii))] += n_;
      incr [ordUxw (quadrant->gt (ii))] += Uxw_;
      incr [ordUyw (quadrant->gt (ii))] += Uyw_;
      incr [ordUxs (quadrant->gt (ii))] += Uxs_;
      incr [ordUys (quadrant->gt (ii))] += Uys_;

      P_plus [ordh   (quadrant->gt (ii))] += std::max(0., h_al  );
      P_plus [ordn   (quadrant->gt (ii))] += std::max(0., n_al  );
      P_plus [ordUxw (quadrant->gt (ii))] += std::max(0., Uxw_al);
      P_plus [ordUyw (quadrant->gt (ii))] += std::max(0., Uyw_al);
      P_plus [ordUxs (quadrant->gt (ii))] += std::max(0., Uxs_al);
      P_plus [ordUys (quadrant->gt (ii))] += std::max(0., Uys_al);

      P_minus [ordh   (quadrant->gt (ii))] += std::min(0., h_al  );
      P_minus [ordn   (quadrant->gt (ii))] += std::min(0., n_al  );
      P_minus [ordUxw (quadrant->gt (ii))] += std::min(0., Uxw_al);
      P_minus [ordUyw (quadrant->gt (ii))] += std::min(0., Uyw_al);
      P_minus [ordUxs (quadrant->gt (ii))] += std::min(0., Uxs_al);
      P_minus [ordUys (quadrant->gt (ii))] += std::min(0., Uys_al);
      
      
    } else {


      incr [ordh   (quadrant->gparent(0,ii))] += h_;
      incr [ordh   (quadrant->gparent(1,ii))] += h_;

      incr [ordn   (quadrant->gparent(0,ii))] += n_;
      incr [ordn   (quadrant->gparent(1,ii))] += n_;
      
      incr [ordUxw (quadrant->gparent(0,ii))] += Uxw_;
      incr [ordUxw (quadrant->gparent(1,ii))] += Uxw_;
      
      incr [ordUyw (quadrant->gparent(0,ii))] += Uyw_;
      incr [ordUyw (quadrant->gparent(1,ii))] += Uyw_;

      incr [ordUxs (quadrant->gparent(0,ii))] += Uxs_;
      incr [ordUxs (quadrant->gparent(1,ii))] += Uxs_;
      
      incr [ordUys (quadrant->gparent(0,ii))] += Uys_;
      incr [ordUys (quadrant->gparent(1,ii))] += Uys_;



      P_plus [ordh   (quadrant->gparent(0,ii))] += std::max(0., h_al  );
      P_plus [ordh   (quadrant->gparent(1,ii))] += std::max(0., h_al  );

      P_plus [ordn   (quadrant->gparent(0,ii))] += std::max(0., n_al  );
      P_plus [ordn   (quadrant->gparent(1,ii))] += std::max(0., n_al  );
      
      P_plus [ordUxw (quadrant->gparent(0,ii))] += std::max(0., Uxw_al);
      P_plus [ordUxw (quadrant->gparent(1,ii))] += std::max(0., Uxw_al);
      
      P_plus [ordUyw (quadrant->gparent(0,ii))] += std::max(0., Uyw_al);
      P_plus [ordUyw (quadrant->gparent(1,ii))] += std::max(0., Uyw_al);

      P_plus [ordUxs (quadrant->gparent(0,ii))] += std::max(0., Uxs_al);
      P_plus [ordUxs (quadrant->gparent(1,ii))] += std::max(0., Uxs_al);
      
      P_plus [ordUys (quadrant->gparent(0,ii))] += std::max(0., Uys_al);
      P_plus [ordUys (quadrant->gparent(1,ii))] += std::max(0., Uys_al);


 
      P_minus [ordh   (quadrant->gparent(0,ii))] += std::min(0., h_al  );
      P_minus [ordh   (quadrant->gparent(1,ii))] += std::min(0., h_al  );

      P_minus [ordn   (quadrant->gparent(0,ii))] += std::min(0., n_al  );
      P_minus [ordn   (quadrant->gparent(1,ii))] += std::min(0., n_al  );
      
      P_minus [ordUxw (quadrant->gparent(0,ii))] += std::min(0., Uxw_al);
      P_minus [ordUxw (quadrant->gparent(1,ii))] += std::min(0., Uxw_al);
      
      P_minus [ordUyw (quadrant->gparent(0,ii))] += std::min(0., Uyw_al);
      P_minus [ordUyw (quadrant->gparent(1,ii))] += std::min(0., Uyw_al);

      P_minus [ordUxs (quadrant->gparent(0,ii))] += std::min(0., Uxs_al);
      P_minus [ordUxs (quadrant->gparent(1,ii))] += std::min(0., Uxs_al);
      
      P_minus [ordUys (quadrant->gparent(0,ii))] += std::min(0., Uys_al);
      P_minus [ordUys (quadrant->gparent(1,ii))] += std::min(0., Uys_al);
      
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

    double hdof_c, ndof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c, P_plus_h_c, P_minus_h_c, P_plus_n_c, P_minus_n_c, P_plus_Uxw_c, P_minus_Uxw_c, P_plus_Uyw_c, P_minus_Uyw_c, P_plus_Uxs_c, P_minus_Uxs_c, P_plus_Uys_c, P_minus_Uys_c;

    if (! quadrant->is_hanging (ii)){
      hdof_c       = sol [ordh     (quadrant->gt (ii))];
      ndof_c       = sol [ordn     (quadrant->gt (ii))];
      Uxwdof_c     = sol [ordUxw   (quadrant->gt (ii))];
      Uywdof_c     = sol [ordUyw   (quadrant->gt (ii))];
      Uxsdof_c     = sol [ordUxs   (quadrant->gt (ii))];
      Uysdof_c     = sol [ordUys   (quadrant->gt (ii))];

      Z_node[ii]   = Z [quadrant->gt (ii)];

      P_plus_h_c   = P_plus [ordh   (quadrant->gt (ii))];
      P_minus_h_c  = P_minus[ordh   (quadrant->gt (ii))];

      P_plus_n_c   = P_plus [ordn   (quadrant->gt (ii))];
      P_minus_n_c  = P_minus[ordn   (quadrant->gt (ii))];

      P_plus_Uxw_c   = P_plus [ordUxw   (quadrant->gt (ii))];
      P_minus_Uxw_c  = P_minus[ordUxw   (quadrant->gt (ii))];

      P_plus_Uyw_c   = P_plus [ordUyw   (quadrant->gt (ii))];
      P_minus_Uyw_c  = P_minus[ordUyw   (quadrant->gt (ii))];

      P_plus_Uxs_c   = P_plus [ordUxs   (quadrant->gt (ii))];
      P_minus_Uxs_c  = P_minus[ordUxs   (quadrant->gt (ii))];

      P_plus_Uys_c   = P_plus [ordUys   (quadrant->gt (ii))];
      P_minus_Uys_c  = P_minus[ordUys   (quadrant->gt (ii))];
      

    } else {
      hdof_c   = .5 * (sol [ordh   (quadrant->gparent(0,ii))] +
                       sol [ordh   (quadrant->gparent(1,ii))]);
      ndof_c   = .5 * (sol [ordn   (quadrant->gparent(0,ii))] +
                       sol [ordn   (quadrant->gparent(1,ii))]);
      Uxwdof_c = .5 * (sol [ordUxw (quadrant->gparent(0,ii))] +
                       sol [ordUxw (quadrant->gparent(1,ii))]);
      Uywdof_c = .5 * (sol [ordUyw (quadrant->gparent(0,ii))] +
                       sol [ordUyw (quadrant->gparent(1,ii))]);
      Uxsdof_c = .5 * (sol [ordUxs (quadrant->gparent(0,ii))] +
                       sol [ordUxs (quadrant->gparent(1,ii))]);
      Uysdof_c = .5 * (sol [ordUys (quadrant->gparent(0,ii))] +
                       sol [ordUys (quadrant->gparent(1,ii))]);

      Z_node[ii] = .5 * (Z [quadrant->gparent(0,ii)] +
                         Z [quadrant->gparent(1,ii)]);

      P_plus_h_c   = .5 * (P_plus [ordh (quadrant->gparent(0,ii))] +
                           P_plus [ordh (quadrant->gparent(1,ii))]);
      P_minus_h_c  = .5 * (P_minus [ordh (quadrant->gparent(0,ii))] +
                           P_minus [ordh (quadrant->gparent(1,ii))]);

      P_plus_n_c   = .5 * (P_plus [ordn (quadrant->gparent(0,ii))] +
                           P_plus [ordn (quadrant->gparent(1,ii))]);
      P_minus_n_c  = .5 * (P_minus [ordn (quadrant->gparent(0,ii))] +
                           P_minus [ordn (quadrant->gparent(1,ii))]);

      P_plus_Uxw_c  = .5 * (P_plus [ordUxw (quadrant->gparent(0,ii))] +
                            P_plus [ordUxw (quadrant->gparent(1,ii))]);
      P_minus_Uxw_c = .5 * (P_minus [ordUxw (quadrant->gparent(0,ii))] +
                            P_minus [ordUxw (quadrant->gparent(1,ii))]);

      P_plus_Uyw_c  = .5 * (P_plus [ordUyw (quadrant->gparent(0,ii))] +
                            P_plus [ordUyw (quadrant->gparent(1,ii))]);
      P_minus_Uyw_c = .5 * (P_minus [ordUyw (quadrant->gparent(0,ii))] +
                            P_minus [ordUyw (quadrant->gparent(1,ii))]);

      P_plus_Uxs_c  = .5 * (P_plus [ordUxs (quadrant->gparent(0,ii))] +
                            P_plus [ordUxs (quadrant->gparent(1,ii))]);
      P_minus_Uxs_c = .5 * (P_minus [ordUxs (quadrant->gparent(0,ii))] +
                            P_minus [ordUxs (quadrant->gparent(1,ii))]);

      P_plus_Uys_c  = .5 * (P_plus [ordUys (quadrant->gparent(0,ii))] +
                            P_plus [ordUys (quadrant->gparent(1,ii))]);
      P_minus_Uys_c = .5 * (P_minus [ordUys (quadrant->gparent(0,ii))] +
                            P_minus [ordUys (quadrant->gparent(1,ii))]);
      
    }

    etadof     [ii] = hdof_c>epsilon ? hdof_c+Z_node[ii] : hdof_c;
    hdof       [ii] = hdof_c;
    ndof       [ii] = ndof_c;
    Uxwdof     [ii] = Uxwdof_c;
    Uywdof     [ii] = Uywdof_c;
    Uxsdof     [ii] = Uxsdof_c;
    Uysdof     [ii] = Uysdof_c;

    P_plus_h_dof [ii] = P_plus_h_c;
    P_minus_h_dof[ii] = P_minus_h_c;

    P_plus_n_dof [ii] = P_plus_n_c;
    P_minus_n_dof[ii] = P_minus_n_c;

    P_plus_Uxw_dof [ii] = P_plus_Uxw_c;
    P_minus_Uxw_dof[ii] = P_minus_Uxw_c;

    P_plus_Uyw_dof [ii] = P_plus_Uyw_c;
    P_minus_Uyw_dof[ii] = P_minus_Uyw_c;

    P_plus_Uxs_dof [ii] = P_plus_Uxs_c;
    P_minus_Uxs_dof[ii] = P_minus_Uxs_c;

    P_plus_Uys_dof [ii] = P_plus_Uys_c;
    P_minus_Uys_dof[ii] = P_minus_Uys_c;

  }
  
  
  
  // compute local extrema
  const auto h_min_cell  = *std::min_element(etadof.begin(), etadof.end());
  const auto h_max_cell  = *std::max_element(etadof.begin(), etadof.end());

  const auto n_min_cell  = *std::min_element(ndof.begin(), ndof.end());
  const auto n_max_cell  = *std::max_element(ndof.begin(), ndof.end());
 
  const auto Uxw_min_cell = *std::min_element(Uxwdof.begin(),  Uxwdof.end() );
  const auto Uxw_max_cell = *std::max_element(Uxwdof.begin(),  Uxwdof.end() );

  const auto Uyw_min_cell = *std::min_element(Uywdof.begin(),  Uywdof.end() );
  const auto Uyw_max_cell = *std::max_element(Uywdof.begin(),  Uywdof.end() );

  const auto Uxs_min_cell = *std::min_element(Uxsdof.begin(),  Uxsdof.end() );
  const auto Uxs_max_cell = *std::max_element(Uxsdof.begin(),  Uxsdof.end() );

  const auto Uys_min_cell = *std::min_element(Uysdof.begin(),  Uysdof.end() );
  const auto Uys_max_cell = *std::max_element(Uysdof.begin(),  Uysdof.end() );  

  bool is_node_in_element = false;

  std::array<double,4> h_min  = {h_min_cell, h_min_cell, h_min_cell, h_min_cell }, h_max  = {h_max_cell, h_max_cell, h_max_cell, h_max_cell },
                       n_min  = {n_min_cell, n_min_cell, n_min_cell, n_min_cell }, n_max  = {n_max_cell, n_max_cell, n_max_cell, n_max_cell },
                       Uxw_min = {Uxw_min_cell,Uxw_min_cell,Uxw_min_cell,Uxw_min_cell}, Uxw_max = {Uxw_max_cell,Uxw_max_cell,Uxw_max_cell,Uxw_max_cell},
                       Uyw_min = {Uyw_min_cell,Uyw_min_cell,Uyw_min_cell,Uyw_min_cell}, Uyw_max = {Uyw_max_cell,Uyw_max_cell,Uyw_max_cell,Uyw_max_cell},
                       Uxs_min = {Uxs_min_cell,Uxs_min_cell,Uxs_min_cell,Uxs_min_cell}, Uxs_max = {Uxs_max_cell,Uxs_max_cell,Uxs_max_cell,Uxs_max_cell},
                       Uys_min = {Uys_min_cell,Uys_min_cell,Uys_min_cell,Uys_min_cell}, Uys_max = {Uys_max_cell,Uys_max_cell,Uys_max_cell,Uys_max_cell};

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

          double h_current_cell, n_current_cell, Uxw_current_cell, Uyw_current_cell, Uxs_current_cell, Uys_current_cell;

          if (! quadrant_nei->is_hanging (jj)){
            h_current_cell   = sol [ordh   (quadrant_nei->gt (jj))] + Z [quadrant->gt (ii)];
            n_current_cell   = sol [ordn   (quadrant_nei->gt (jj))];
            Uxw_current_cell = sol [ordUxw (quadrant_nei->gt (jj))];
            Uyw_current_cell = sol [ordUyw (quadrant_nei->gt (jj))];
            Uxs_current_cell = sol [ordUxs (quadrant_nei->gt (jj))];
            Uys_current_cell = sol [ordUys (quadrant_nei->gt (jj))];
          } else {
            h_current_cell  = .5 * (sol [ordh  (quadrant_nei->gparent(0,jj))] +
                                    sol [ordh  (quadrant_nei->gparent(1,jj))]) + 
                              .5 * (Z [quadrant->gparent(0,ii)] +
                                    Z [quadrant->gparent(1,ii)]);
            n_current_cell  = .5 * (sol [ordn  (quadrant_nei->gparent(0,jj))] +
                                    sol [ordn  (quadrant_nei->gparent(1,jj))]);
            Uxw_current_cell = .5 * (sol [ordUxw (quadrant_nei->gparent(0,jj))] +
                                    sol [ordUxw (quadrant_nei->gparent(1,jj))]);
            Uyw_current_cell = .5 * (sol [ordUyw (quadrant_nei->gparent(0,jj))] +
                                    sol [ordUyw (quadrant_nei->gparent(1,jj))]);
            Uxs_current_cell = .5 * (sol [ordUxs (quadrant_nei->gparent(0,jj))] +
                                    sol [ordUxs (quadrant_nei->gparent(1,jj))]);
            Uys_current_cell = .5 * (sol [ordUys (quadrant_nei->gparent(0,jj))] +
                                    sol [ordUys (quadrant_nei->gparent(1,jj))]);
          }

          h_min[ii]  = std::min(h_min[ii],  h_current_cell);
          h_max[ii]  = std::max(h_max[ii],  h_current_cell);

          n_min[ii]  = std::min(n_min[ii],  n_current_cell);
          n_max[ii]  = std::max(n_max[ii],  n_current_cell);

          Uxw_min[ii] = std::min(Uxw_min[ii], Uxw_current_cell);
          Uxw_max[ii] = std::max(Uxw_max[ii], Uxw_current_cell);

          Uyw_min[ii] = std::min(Uyw_min[ii], Uyw_current_cell);
          Uyw_max[ii] = std::max(Uyw_max[ii], Uyw_current_cell);

          Uxs_min[ii] = std::min(Uxs_min[ii], Uxs_current_cell);
          Uxs_max[ii] = std::max(Uxs_max[ii], Uxs_current_cell);

          Uys_min[ii] = std::min(Uys_min[ii], Uys_current_cell);
          Uys_max[ii] = std::max(Uys_max[ii], Uys_current_cell);

        }

        is_node_in_element = false;
      }


    }
  }



  
  // compute flux correction
  double phi_cell_h = 1., phi_cell_n = 1., phi_cell_Uxw = 1., phi_cell_Uyw = 1., phi_cell_Uxs = 1., phi_cell_Uys = 1.;
  for (int ii = 0; ii < 4; ++ii){

    const auto & flux_on_the_node_h   = incr_anti_diff[ordh  (index_quadrant)][ii];
    const auto & flux_on_the_node_n   = incr_anti_diff[ordn  (index_quadrant)][ii];
    const auto & flux_on_the_node_Uxw = incr_anti_diff[ordUxw(index_quadrant)][ii];
    const auto & flux_on_the_node_Uyw = incr_anti_diff[ordUyw(index_quadrant)][ii];
    const auto & flux_on_the_node_Uxs = incr_anti_diff[ordUxs(index_quadrant)][ii];
    const auto & flux_on_the_node_Uys = incr_anti_diff[ordUys(index_quadrant)][ii];


    const auto speed = max_eigen(hdof[ii], ndof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii]);
    const auto vel_rusanov_cell_x = speed[0]; //hpoint>epsilon ? (std::abs(Uxdof[ii]/hpoint)+celerity) : 0.;
    const auto vel_rusanov_cell_y = speed[1]; //hpoint>epsilon ? (std::abs(Uydof[ii]/hpoint)+celerity) : 0.;

    const auto vel_square_rusanov_cell = vel_rusanov_cell_x * vel_rusanov_cell_y;

    flux_limiter(h_min  [ii], h_max  [ii], etadof [ii], P_plus_h_dof  [ii], P_minus_h_dof   [ii], flux_on_the_node_h,   vel_square_rusanov_cell, phi_cell_h  );
    flux_limiter(n_min  [ii], n_max  [ii], ndof   [ii], P_plus_n_dof  [ii], P_minus_n_dof   [ii], flux_on_the_node_n,   vel_square_rusanov_cell, phi_cell_n  );
    flux_limiter(Uxw_min[ii], Uxw_max[ii], Uxwdof [ii], P_plus_Uxw_dof[ii], P_minus_Uxw_dof [ii], flux_on_the_node_Uxw, vel_square_rusanov_cell, phi_cell_Uxw);
    flux_limiter(Uyw_min[ii], Uyw_max[ii], Uywdof [ii], P_plus_Uyw_dof[ii], P_minus_Uyw_dof [ii], flux_on_the_node_Uyw, vel_square_rusanov_cell, phi_cell_Uyw);
    flux_limiter(Uxs_min[ii], Uxs_max[ii], Uxsdof [ii], P_plus_Uxs_dof[ii], P_minus_Uxs_dof [ii], flux_on_the_node_Uxs, vel_square_rusanov_cell, phi_cell_Uxs);
    flux_limiter(Uys_min[ii], Uys_max[ii], Uysdof [ii], P_plus_Uys_dof[ii], P_minus_Uys_dof [ii], flux_on_the_node_Uys, vel_square_rusanov_cell, phi_cell_Uys);
  }

  //std::cout << phi_cell_h << " " << phi_cell_Ux << " " << phi_cell_Uy << std::endl;

  //<phi_cell_h = 1., phi_cell_Ux = 1., phi_cell_Uy = 1.; 
  //std::cout << "we put limiter equal to one" << std::endl;



  for (int ii = 0; ii < 4; ++ii){

    //std::cout << phi_cell_h << " " << phi_cell_Ux << " " << phi_cell_Uy << std::endl;

    const auto flux_on_the_node_h   = incr_anti_diff[ordh  (index_quadrant)][ii]*phi_cell_h;
    const auto flux_on_the_node_n   = incr_anti_diff[ordn  (index_quadrant)][ii]*phi_cell_n;
    const auto flux_on_the_node_Uxw = incr_anti_diff[ordUxw(index_quadrant)][ii]*phi_cell_Uxw;
    const auto flux_on_the_node_Uyw = incr_anti_diff[ordUyw(index_quadrant)][ii]*phi_cell_Uyw;
    const auto flux_on_the_node_Uxs = incr_anti_diff[ordUxs(index_quadrant)][ii]*phi_cell_Uxs;
    const auto flux_on_the_node_Uys = incr_anti_diff[ordUys(index_quadrant)][ii]*phi_cell_Uys;

    if (! quadrant->is_hanging (ii)){

      incr [ordh   (quadrant->gt (ii))] += flux_on_the_node_h;
      incr [ordn   (quadrant->gt (ii))] += flux_on_the_node_n;
      incr [ordUxw (quadrant->gt (ii))] += flux_on_the_node_Uxw;
      incr [ordUyw (quadrant->gt (ii))] += flux_on_the_node_Uyw;
      incr [ordUxs (quadrant->gt (ii))] += flux_on_the_node_Uxs;
      incr [ordUys (quadrant->gt (ii))] += flux_on_the_node_Uys;

    } else {

      incr [ordh   (quadrant->gparent(0,ii))] += flux_on_the_node_h; 
      incr [ordh   (quadrant->gparent(1,ii))] += flux_on_the_node_h;

      incr [ordn   (quadrant->gparent(0,ii))] += flux_on_the_node_n; 
      incr [ordn   (quadrant->gparent(1,ii))] += flux_on_the_node_n;
      
      incr [ordUxw (quadrant->gparent(0,ii))] += flux_on_the_node_Uxw;
      incr [ordUxw (quadrant->gparent(1,ii))] += flux_on_the_node_Uxw;
      
      incr [ordUyw (quadrant->gparent(0,ii))] += flux_on_the_node_Uyw;
      incr [ordUyw (quadrant->gparent(1,ii))] += flux_on_the_node_Uyw;

      incr [ordUxs (quadrant->gparent(0,ii))] += flux_on_the_node_Uxs;
      incr [ordUxs (quadrant->gparent(1,ii))] += flux_on_the_node_Uxs;
      
      incr [ordUys (quadrant->gparent(0,ii))] += flux_on_the_node_Uys;
      incr [ordUys (quadrant->gparent(1,ii))] += flux_on_the_node_Uys;

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
TG2_scheme::h_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys) 
{ 
  // A flux-limiting wetting–drying method for finite-element shallow-water models, with application to the Scheldt Estuary
  return (h>epsilon ? Uxw+Uxs : 0.); 
}

double
TG2_scheme::h_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  // A flux-limiting wetting–drying method for finite-element shallow-water models, with application to the Scheldt Estuary
  return (h>epsilon ? Uyw+Uys : 0.); 
}

double
TG2_scheme::n_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  return (0.); 
}

double
TG2_scheme::n_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  return (0.); 
}

double
TG2_scheme::Uwx_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto vel_x = (n*h)>epsilon ? Uxw/(n*h) : 0.;
  return (Uxw*vel_x + grav*h*h/2.*n); 
}
 
double
TG2_scheme::Uwx_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return ((n*h)>epsilon ? Uyw*Uxw/(n*h) : 0.); }

double
TG2_scheme::Uwy_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return ((n*h)>epsilon ? Uyw*Uxw/(n*h) : 0.); }

double
TG2_scheme::Uwy_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto vel_y = (n*h)>epsilon ? Uyw/(n*h) : 0.;
  return (Uyw*vel_y + grav*h*h/2.*n); 
}


double
TG2_scheme::Usx_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto vel_x = ((1.-n)*h)>epsilon ? Uxs/((1.-n)*h) : 0.;
  return (Uxs*vel_x + grav*h*h/2.*(1.-n)); 
}
 
double
TG2_scheme::Usx_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (((1.-n)*h)>epsilon ? Uys*Uxs/((1.-n)*h) : 0.); }

double
TG2_scheme::Usy_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (((1.-n)*h)>epsilon ? Uys*Uxs/((1.-n)*h) : 0.); }

double
TG2_scheme::Usy_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto vel_y = ((1.-n)*h)>epsilon ? Uys/((1.-n)*h) : 0.;
  return (Uys*vel_y + grav*h*h/2.*(1.-n)); 
}




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
TG2_scheme::h_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto Ux = Uxw + Uxs; 
  const auto Uy = Uyw + Uys;
  return (erosion_coefficient*std::sqrt(Ux*Ux+Uy*Uy)); 
}

double
TG2_scheme::Uxw_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const double hw = n*h;
  const double hs = (1.-n)*h;

  const double vel_w_x = hw>epsilon ? Uxw/hw : 0.;
  const double vel_s_x = hs>epsilon ? Uxs/hs : 0.;

  const double C_d = hw>epsilon && h>epsilon ? n*(1.-n)/std::pow(n, m_coeff)/terminal_velocity*(density_s/density_w-1.)*grav : 0.;
  const double R_x = C_d*(vel_w_x-vel_s_x);

  return ( -h*R_x );
}

double
TG2_scheme::Uyw_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const double hw = n*h;
  const double hs = (1.-n)*h;

  const double vel_w_y = hw>epsilon ? Uyw/hw : 0.;
  const double vel_s_y = hs>epsilon ? Uys/hs : 0.;

  const double C_d = hw>epsilon && h>epsilon ? n*(1.-n)/std::pow(n, m_coeff)/terminal_velocity*(density_s/density_w-1.)*grav : 0.;
  const double R_y = C_d*(vel_w_y-vel_s_y);

  return ( -h*R_y );
}

double
TG2_scheme::Uxs_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const double bed_pressure = grav*h*(1.-n)*(1.-density_w/density_s); 
  const auto Ux = Uxw + Uxs; 
  const auto Uy = Uyw + Uys;
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::abs( vel_x );

  const double hw = n*h;
  const double hs = (1.-n)*h;

  const double vel_w_x = hw>epsilon ? Uxw/hw : 0.;
  const double vel_s_x = hs>epsilon ? Uxs/hs : 0.;

  const double density = (1.-n) + n*density_w/density_s;

  const double C_d = hw>epsilon && h>epsilon ? n*(1.-n)/std::pow(n, m_coeff)/terminal_velocity*(1.-density_w/density_s)*grav : 0.;
  const double R_x = C_d*(vel_w_x-vel_s_x);

  //std::cout << h << " " << dhdx << " " << dZdx << " " << grav*h*(dZdx+dhdx) << std::endl;

  const double vel_x_sign = std::abs(Ux)>tolerance_sign ? Ux/std::abs(Ux) : Ux/tolerance_sign;
  //const double vel_x_sign = abs_vel>tolerance_sign ? vel_x/abs_vel : 0.;
  //const double vel_x_sign = abs_vel>tolerance_sign ? vel_x/abs_vel : vel_x/tolerance_sign;

  //const double bed_fric_contr = is_bed_friction ? vel_x_sign*(grav*abs_vel*abs_vel/turbulence_coeff + bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;

  const double bed_fric_contr_one = h*h>epsilon ? density*Ux*grav*std::abs(Ux)/turbulence_coeff/h/h : 0.; //is_bed_friction ? vel_x*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = vel_x_sign*bed_pressure*std::tan(bed_friction_angle_rad);

  return ( - bed_fric_contr_one - bed_fric_contr_two + h*R_x );
}

double
TG2_scheme::Uys_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const double bed_pressure = grav*h*(1.-n)*(1.-density_w/density_s);
  const auto Ux = Uxw + Uxs; 
  const auto Uy = Uyw + Uys;
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::abs( vel_y );

  const double hw = n*h;
  const double hs = (1.-n)*h;

  const double vel_w_y = hw>epsilon ? Uyw/hw : 0.;
  const double vel_s_y = hs>epsilon ? Uys/hs : 0.;

  const double density = (1.-n) + n*density_w/density_s;

  const double C_d = hw>epsilon && h>epsilon ? n*(1.-n)/std::pow(n, m_coeff)/terminal_velocity*(1.-density_w/density_s)*grav : 0.;

  const double R_y = C_d*(vel_w_y-vel_s_y);

  const double vel_y_sign = std::abs(Uy)>tolerance_sign ? Uy/std::abs(Uy) : Uy/tolerance_sign;
  //const double vel_y_sign = abs_vel>tolerance_sign ? vel_y/abs_vel : vel_y/tolerance_sign;
  //const double vel_y_sign = (vel_y > ) ? 1.0 : (vel_y < 0) ? -1.0 : 0.0;

  //const double bed_fric_contr = is_bed_friction ? vel_y_sign*(grav*abs_vel*abs_vel/turbulence_coeff + bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;

  const double bed_fric_contr_one = h*h>epsilon ? density*Uy*grav*std::abs(Uy)/turbulence_coeff/h/h : 0.; //is_bed_friction ? vel_y*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = vel_y_sign*bed_pressure*std::tan(bed_friction_angle_rad);

  //std::cout << dZdy << std::endl; 

  return ( - bed_fric_contr_one - bed_fric_contr_two + h*R_y );
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


double
TG2_scheme::signum (const double& x)
{ return ((x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0); }


