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
             Q1& P_plus,
             Q1& P_minus,
             Q0& sol_onehalf,
             Q1& mass,
             const ordering& ohw,
             const ordering& ohs,
             const ordering& oUxw,
             const ordering& oUyw,
             const ordering& oUxs,
             const ordering& oUys,
             Q1& Z,
             Q0& Z_onehalf,
             const double& DELTAT,
             const double& h_min,
             const bool& is_non_reflBC,
             const bool& is_bed_friction,
             const double& grav,
             const double& density_w,
             const double& density_s,
             const double& turbulence_coeff,
             const double& bed_friction_angle_rad,
             const double& erosion_coefficient,
             const double& m_coeff,
             const double& terminal_velocity,
             std::vector<double>& slope_x,
             std::vector<double>& slope_y)
: sol(sol), 
sold(sold), 
soldd(soldd), 
sold_rkc(sold_rkc), 
soldd_rkc(soldd_rkc), 
sol_ini_rkc(sol_ini_rkc), 
incr(incr), 
incr_initial_source(incr_initial_source),
incr_source(incr_source),
incr_anti_diff(incr_anti_diff), 
P_plus(P_plus), 
P_minus(P_minus), 
sol_onehalf(sol_onehalf), 
mass(mass), 
ordhw(ohw), 
ordhs(ohs), 
ordUxw(oUxw), 
ordUyw(oUyw), 
ordUxs(oUxs), 
ordUys(oUys), 
Z(Z), 
Z_onehalf(Z_onehalf), 
DELTAT(DELTAT), 
epsilon(h_min), 
is_non_reflBC(is_non_reflBC), 
is_bed_friction(is_bed_friction), 
grav(grav),
density_w(density_w), 
density_s(density_s), 
turbulence_coeff(turbulence_coeff), 
bed_friction_angle_rad(bed_friction_angle_rad), 
erosion_coefficient(erosion_coefficient),
m_coeff(m_coeff), 
terminal_velocity(terminal_velocity),
slope_x(slope_x),
slope_y(slope_y)
{ }
 

std::array<double,2>
TG2_scheme::max_eigen (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  double h = hw+hs,

         vel_x  = h >epsilon ? (Uxw+Uxs)/h : 0., 
         vel_y  = h >epsilon ? (Uyw+Uys)/h : 0., 
         velw_x = hw>epsilon ? Uxw/hw : 0., 
         velw_y = hw>epsilon ? Uyw/hw : 0., 
         vels_x = hs>epsilon ? Uxs/hs : 0., 
         vels_y = hs>epsilon ? Uys/hs : 0., 

         celerity = std::sqrt(grav*h);


  const auto lambda_x = h>epsilon ? std::max(std::abs(velw_x), std::abs(vels_x))+celerity : 0.;
  const auto lambda_y = h>epsilon ? std::max(std::abs(velw_y), std::abs(vels_y))+celerity : 0.;

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
      hwdof[ii]  = sol [ordhw  (quadrant->gt (ii) )];
      hsdof[ii]  = sol [ordhs  (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 
    }
    else
    {
      hwdof[ii]  = .5 * (sol [ordhw  (quadrant->gparent (0, ii) )] +
                         sol [ordhw  (quadrant->gparent (1, ii) )]);
      hsdof[ii]  = .5 * (sol [ordhs  (quadrant->gparent (0, ii) )] +
                         sol [ordhs  (quadrant->gparent (1, ii) )]);
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

    const auto lambdas = max_eigen (hwdof[ii], hsdof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii]);
    const auto & vel_rusanov_cell_x = lambdas[0];
    const auto & vel_rusanov_cell_y = lambdas[1];

    const auto h_c = hwdof[ii]+hsdof[ii];

    const auto dtoptx = h_c>epsilon ? Dx/vel_rusanov_cell_x : DELTAT;
    const auto dtopty = h_c>epsilon ? Dy/vel_rusanov_cell_y : DELTAT;
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
    double hwdof, hwdof_old, hwdof_oldd, hsdof, hsdof_old, hsdof_oldd, hdof, hdof_old, hdof_oldd, dh_t, h1, h2, h3, a_coeff, b_coeff;
    if (! quadrant->is_hanging (ii) )
    {
      hwdof      = sol   [ordhw  (quadrant->gt (ii) )];
      hwdof_old  = sold  [ordhw  (quadrant->gt (ii) )];
      hwdof_oldd = soldd [ordhw  (quadrant->gt (ii) )];

      hsdof      = sol   [ordhs  (quadrant->gt (ii) )];
      hsdof_old  = sold  [ordhs  (quadrant->gt (ii) )];
      hsdof_oldd = soldd [ordhs  (quadrant->gt (ii) )];
    }
    else
    {
      hwdof      = .5 * (sol   [ordhw  (quadrant->gparent (0, ii) )] +
                         sol   [ordhw  (quadrant->gparent (1, ii) )]);
      hwdof_old  = .5 * (sold  [ordhw  (quadrant->gparent (0, ii) )] +
                         sold  [ordhw  (quadrant->gparent (1, ii) )]);
      hwdof_oldd = .5 * (soldd [ordhw  (quadrant->gparent (0, ii) )] +
                         soldd [ordhw  (quadrant->gparent (1, ii) )]);

      hsdof      = .5 * (sol   [ordhs  (quadrant->gparent (0, ii) )] +
                         sol   [ordhs  (quadrant->gparent (1, ii) )]);
      hsdof_old  = .5 * (sold  [ordhs  (quadrant->gparent (0, ii) )] +
                         sold  [ordhs  (quadrant->gparent (1, ii) )]);
      hsdof_oldd = .5 * (soldd [ordhs  (quadrant->gparent (0, ii) )] +
                         soldd [ordhs  (quadrant->gparent (1, ii) )]);
    }
    hdof = hwdof+hsdof;
    hdof_old = hwdof_old+hsdof_old;
    hdof_oldd = hwdof_oldd+hsdof_oldd;

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
  
  const auto & index_quadrant = quadrant->get_global_quad_idx (); 
  
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii);
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;


  double hw_cell_average = 0., hs_cell_average = 0., Uxw_cell_average = 0., Uyw_cell_average = 0., Uxs_cell_average = 0., Uys_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hwdof [ii] = sol [ordhw  (quadrant->gt (ii) )];
      hsdof [ii] = sol [ordhs  (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 

      Z_node[ii]  = Z [quadrant->gt (ii)];
    }
    else
    {
      hwdof [ii] = .5 * (sol [ordhw  (quadrant->gparent (0, ii) )] +
                         sol [ordhw  (quadrant->gparent (1, ii) )]);
      hsdof [ii] = .5 * (sol [ordhs  (quadrant->gparent (0, ii) )] +
                         sol [ordhs  (quadrant->gparent (1, ii) )]);
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

    const auto & hwdof_c  = hwdof [ii];
    const auto & hsdof_c  = hsdof [ii];
    const auto & Uxwdof_c = Uxwdof[ii];
    const auto & Uywdof_c = Uywdof[ii];
    const auto & Uxsdof_c = Uxsdof[ii];
    const auto & Uysdof_c = Uysdof[ii];

    hw_cell_average  += hwdof_c;
    hs_cell_average  += hsdof_c;
    Uxw_cell_average += Uxwdof_c;
    Uyw_cell_average += Uywdof_c;
    Uxs_cell_average += Uxsdof_c;
    Uys_cell_average += Uysdof_c;
    
    
    fluxx_hw_node [ii] = hw_flux_formula_x   (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_hw_node [ii] = hw_flux_formula_y   (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_hs_node [ii] = hs_flux_formula_x   (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_hs_node [ii] = hs_flux_formula_y   (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Uxw_node[ii] = Uxw_flux_formula_x  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Uxw_node[ii] = Uxw_flux_formula_y  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Uyw_node[ii] = Uyw_flux_formula_x  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Uyw_node[ii] = Uyw_flux_formula_y  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Uxs_node[ii] = Uxs_flux_formula_x  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Uxs_node[ii] = Uxs_flux_formula_y  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxx_Uys_node[ii] = Uys_flux_formula_x  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
    fluxy_Uys_node[ii] = Uys_flux_formula_y  (hwdof_c, hsdof_c, Uxwdof_c, Uywdof_c, Uxsdof_c, Uysdof_c);
  }
  hw_cell_average  /= 4.;
  hs_cell_average  /= 4.;
  Uxw_cell_average /= 4.;
  Uyw_cell_average /= 4.;
  Uxs_cell_average /= 4.;
  Uys_cell_average /= 4.;
  

  
  const auto div_Fhw_x = .5*((fluxx_hw_node[1]-fluxx_hw_node[0]) + (fluxx_hw_node[3]-fluxx_hw_node[2]));
  const auto div_Fhw_y = .5*((fluxy_hw_node[2]-fluxy_hw_node[0]) + (fluxy_hw_node[3]-fluxy_hw_node[1]));
  const auto div_Fhw_cell = Dy*div_Fhw_x + Dx*div_Fhw_y;

  const auto div_Fhs_x = .5*((fluxx_hs_node[1]-fluxx_hs_node[0]) + (fluxx_hs_node[3]-fluxx_hs_node[2]));
  const auto div_Fhs_y = .5*((fluxy_hs_node[2]-fluxy_hs_node[0]) + (fluxy_hs_node[3]-fluxy_hs_node[1]));
  const auto div_Fhs_cell = Dy*div_Fhs_x + Dx*div_Fhs_y;
  
  const auto div_FUxw_x = .5*((fluxx_Uxw_node[1]-fluxx_Uxw_node[0]) + (fluxx_Uxw_node[3]-fluxx_Uxw_node[2]));
  const auto div_FUxw_y = .5*((fluxy_Uxw_node[2]-fluxy_Uxw_node[0]) + (fluxy_Uxw_node[3]-fluxy_Uxw_node[1]));
  const auto div_FUxw_cell = Dy*div_FUxw_x + Dx*div_FUxw_y;
  
  const auto div_FUyw_x = .5*((fluxx_Uyw_node[1]-fluxx_Uyw_node[0]) + (fluxx_Uyw_node[3]-fluxx_Uyw_node[2]));
  const auto div_FUyw_y = .5*((fluxy_Uyw_node[2]-fluxy_Uyw_node[0]) + (fluxy_Uyw_node[3]-fluxy_Uyw_node[1]));
  const auto div_FUyw_cell = Dy*div_FUyw_x + Dx*div_FUyw_y;

  const auto div_FUxs_x = .5*((fluxx_Uxs_node[1]-fluxx_Uxs_node[0]) + (fluxx_Uxs_node[3]-fluxx_Uxs_node[2]));
  const auto div_FUxs_y = .5*((fluxy_Uxs_node[2]-fluxy_Uxs_node[0]) + (fluxy_Uxs_node[3]-fluxy_Uxs_node[1]));
  const auto div_FUxs_cell = Dy*div_FUxs_x + Dx*div_FUxs_y;
  
  const auto div_FUys_x = .5*((fluxx_Uys_node[1]-fluxx_Uys_node[0]) + (fluxx_Uys_node[3]-fluxx_Uys_node[2]));
  const auto div_FUys_y = .5*((fluxy_Uys_node[2]-fluxy_Uys_node[0]) + (fluxy_Uys_node[3]-fluxy_Uys_node[1]));
  const auto div_FUys_cell = Dy*div_FUys_x + Dx*div_FUys_y;


  const auto slope_x_c = ((Z_node[1] - Z_node[0]) + (Z_node[3] - Z_node[2]))/Dx/2.;
  const auto slope_y_c = ((Z_node[2] - Z_node[0]) + (Z_node[3] - Z_node[1]))/Dy/2.;

  const auto grad_hw_x = ((hwdof[1] - hwdof[0]) + (hwdof[3] - hwdof[2]))/Dx/2.;
  const auto grad_hw_y = ((hwdof[2] - hwdof[0]) + (hwdof[3] - hwdof[1]))/Dy/2.;

  const auto grad_hs_x = ((hsdof[1] - hsdof[0]) + (hsdof[3] - hsdof[2]))/Dx/2.;
  const auto grad_hs_y = ((hsdof[2] - hsdof[0]) + (hsdof[3] - hsdof[1]))/Dy/2.;


  const double tau =(dt + dt_old)*.5*.5;


  Z_onehalf[index_quadrant] = (Z_node[0]+Z_node[1]+Z_node[2]+Z_node[3])*.25;

  auto & hw_c  = sol_onehalf[ordhw   (index_quadrant)];
  auto & hs_c  = sol_onehalf[ordhs   (index_quadrant)];
  auto & Uxw_c = sol_onehalf[ordUxw  (index_quadrant)];
  auto & Uyw_c = sol_onehalf[ordUyw  (index_quadrant)];
  auto & Uxs_c = sol_onehalf[ordUxs  (index_quadrant)];
  auto & Uys_c = sol_onehalf[ordUys  (index_quadrant)];

  const auto v_hw   = hw_cell_average  - tau *  div_Fhw_cell /area;
  const auto v_hs   = hs_cell_average  - tau *  div_Fhs_cell /area;
  const auto v_Ux_w = Uxw_cell_average - tau * (div_FUxw_cell/area - src_slope_formula (hw_cell_average, slope_x_c) +         grav*hw_cell_average*grad_hs_x);
  const auto v_Uy_w = Uyw_cell_average - tau * (div_FUyw_cell/area - src_slope_formula (hw_cell_average, slope_y_c) +         grav*hw_cell_average*grad_hs_y);
  const auto v_Ux_s = Uxs_cell_average - tau * (div_FUxs_cell/area - src_slope_formula (hs_cell_average, slope_x_c) + r_coeff*grav*hs_cell_average*grad_hw_x);
  const auto v_Uy_s = Uys_cell_average - tau * (div_FUys_cell/area - src_slope_formula (hs_cell_average, slope_y_c) + r_coeff*grav*hs_cell_average*grad_hw_y);

  hw_c  = v_hw;
  hs_c  = v_hs;
  Uxw_c = v_Ux_w;
  Uyw_c = v_Uy_w;
  Uxs_c = v_Ux_s;
  Uys_c = v_Uy_s;

  hw_c = hw_c>0 ? hw_c : 0.;
  hs_c = hs_c>0 ? hs_c : 0.;

  const double h_c = hw_c+hs_c;

  const double n_c  = h_c>epsilon ? hw_c/h_c : 0.;
  const double ns_c = h_c>epsilon ? hs_c/h_c : 0.;

  const double density = ns_c*density_s + n_c*density_w;
  const double density_prime = ns_c*(density_s-density_w);


  // solve non-linearities
  const double tolerance = 1.e-4;
  const int Nmax = 1.e3;
  int count = -1;
  double error = tolerance + 1.; 
  while (count++<Nmax && error>tolerance)
  {

    const auto U_tot_x = Uxs_c + Uxw_c;
    const auto U_tot_y = Uys_c + Uyw_c;

    const double Cw_d = (hw_c>epsilon && h_c>epsilon) ? ns_c/terminal_velocity/std::pow(n_c, m_coeff)*(density_s-density_w)*grav : 0.;
    const double Cs_d = (hw_c>epsilon && h_c>epsilon) ? n_c /terminal_velocity/std::pow(n_c, m_coeff)*(density_s-density_w)*grav : 0.;

    const double delta_x = std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign;
    const double delta_y = std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign;

    const double fric_x = (h_c*h_c)>epsilon ? (density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_x)) : 0.;
    const double fric_y = (h_c*h_c)>epsilon ? (density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_y)) : 0.;

    //if (Cw_d!=0)
    //std::cout << mu_tilde_vect[1]*dt* Cw_d/density_w << std::endl;

    const double big_A  = 1. + tau* Cw_d/density_w;
    const double big_B  =    - tau* Cs_d/density_w;
    const double big_Cx =      tau* (is_bed_friction ? (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*delta_x + fric_x) : 0.) /density_s - tau*Cw_d/density_s;
    const double big_Cy =      tau* (is_bed_friction ? (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*delta_y + fric_y) : 0.) /density_s - tau*Cw_d/density_s;
    const double big_Dx = 1. + tau* (is_bed_friction ? (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*delta_x + fric_x) : 0.) /density_s + tau*Cs_d/density_s;
    const double big_Dy = 1. + tau* (is_bed_friction ? (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*delta_y + fric_y) : 0.) /density_s + tau*Cs_d/density_s;


    const double big_detx = big_A*big_Dx-big_B*big_Cx;
    const double big_dety = big_A*big_Dy-big_B*big_Cy;

    const double f_Uwx = - Uxw_c + v_Ux_w + tau*Uxw_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
    const double f_Uwy = - Uyw_c + v_Uy_w + tau*Uyw_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
    const double f_Usx = - Uxs_c + v_Ux_s + tau*Uxs_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
    const double f_Usy = - Uys_c + v_Uy_s + tau*Uys_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);

    const double delta_Uwx = (f_Uwx*big_Dx-f_Usx*big_B )/big_detx;
    const double delta_Usx = (f_Usx*big_A -f_Uwx*big_Cx)/big_detx;

    const double delta_Uwy = (f_Uwy*big_Dy-f_Usy*big_B )/big_dety;
    const double delta_Usy = (f_Usy*big_A -f_Uwy*big_Cy)/big_dety;

    error = std::sqrt(delta_Uwx*delta_Uwx + delta_Uwy*delta_Uwy + delta_Usx*delta_Usx + delta_Usy*delta_Usy);

    Uxw_c += delta_Uwx;
    Uyw_c += delta_Uwy;
    Uxs_c += delta_Usx;
    Uys_c += delta_Usy;

    //std::cout << count << " " << error << std::endl;

  }

}


void
TG2_scheme::solve_non_lin(const int& kk)
{
  auto & hw_c  = sol.get_owned_data ()[kk  ];
  auto & hs_c  = sol.get_owned_data ()[kk+1];
  auto & Uxw_c = sol.get_owned_data ()[kk+2];
  auto & Uyw_c = sol.get_owned_data ()[kk+3];
  auto & Uxs_c = sol.get_owned_data ()[kk+4];
  auto & Uys_c = sol.get_owned_data ()[kk+5];

  double tau = (dt+dt_old)*.5*.5;

  // solve non-linearities like the first step of the TG2 method to get the complete low order solution,
  // terminated this part add the nodal correction to get the hyperbolicity,
  const auto v_hw   = hw_c  + 2.*tau*incr.get_owned_data ()[kk  ]/mass.get_owned_data ()[kk  ];
  const auto v_hs   = hs_c  + 2.*tau*incr.get_owned_data ()[kk+1]/mass.get_owned_data ()[kk+1];
  const auto v_Ux_w = Uxw_c + 2.*tau*incr.get_owned_data ()[kk+2]/mass.get_owned_data ()[kk+2] + tau*Uxw_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
  const auto v_Uy_w = Uyw_c + 2.*tau*incr.get_owned_data ()[kk+3]/mass.get_owned_data ()[kk+3] + tau*Uyw_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
  const auto v_Ux_s = Uxs_c + 2.*tau*incr.get_owned_data ()[kk+4]/mass.get_owned_data ()[kk+4] + tau*Uxs_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
  const auto v_Uy_s = Uys_c + 2.*tau*incr.get_owned_data ()[kk+5]/mass.get_owned_data ()[kk+5] + tau*Uys_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);

  hw_c  = v_hw;
  hs_c  = v_hs;
  Uxw_c = v_Ux_w;
  Uyw_c = v_Uy_w;
  Uxs_c = v_Ux_s;
  Uys_c = v_Uy_s;

  hw_c = hw_c>0 ? hw_c : 0.;
  hs_c = hs_c>0 ? hs_c : 0.;

  const double h_c = hw_c+hs_c;

  const double n_c  = h_c>epsilon ? hw_c/h_c : 0.;
  const double ns_c = h_c>epsilon ? hs_c/h_c : 0.;

  const double density = ns_c*density_s + n_c*density_w;
  const double density_prime = ns_c*(density_s-density_w);


  // solve non-linearities
  const double tolerance = 1.e-4;
  const int Nmax = 1.e3;
  int count = -1;
  double error = tolerance + 1.; 
  while (count++<Nmax && error>tolerance)
  {
    const auto U_tot_x = Uxs_c + Uxw_c;
    const auto U_tot_y = Uys_c + Uyw_c;

    const double Cw_d = (hw_c>epsilon && h_c>epsilon) ? ns_c/terminal_velocity/std::pow(n_c, m_coeff)*(density_s-density_w)*grav : 0.;
    const double Cs_d = (hw_c>epsilon && h_c>epsilon) ? n_c /terminal_velocity/std::pow(n_c, m_coeff)*(density_s-density_w)*grav : 0.;

    const double delta_x = std::abs(U_tot_x)>tolerance_sign ? 0. : 1./tolerance_sign;
    const double delta_y = std::abs(U_tot_y)>tolerance_sign ? 0. : 1./tolerance_sign;

    const double fric_x = (h_c*h_c)>epsilon ? (density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_x)) : 0.;
    const double fric_y = (h_c*h_c)>epsilon ? (density*grav/turbulence_coeff/(h_c*h_c)*2.*std::abs(U_tot_y)) : 0.;

    //if (Cw_d!=0)
    //std::cout << mu_tilde_vect[1]*dt* Cw_d/density_w << std::endl;

    const double big_A  = 1. + tau* Cw_d/density_w;
    const double big_B  =    - tau* Cs_d/density_w;
    const double big_Cx =      tau* (is_bed_friction ? (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*delta_x + fric_x) : 0.) /density_s - tau*Cw_d/density_s;
    const double big_Cy =      tau* (is_bed_friction ? (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*delta_y + fric_y) : 0.) /density_s - tau*Cw_d/density_s;
    const double big_Dx = 1. + tau* (is_bed_friction ? (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*delta_x + fric_x) : 0.) /density_s + tau*Cs_d/density_s;
    const double big_Dy = 1. + tau* (is_bed_friction ? (density_prime*grav*h_c*std::tan(bed_friction_angle_rad)*delta_y + fric_y) : 0.) /density_s + tau*Cs_d/density_s;


    const double big_detx = big_A*big_Dx-big_B*big_Cx;
    const double big_dety = big_A*big_Dy-big_B*big_Cy;

    const double f_Uwx = - Uxw_c + v_Ux_w + tau*Uxw_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
    const double f_Uwy = - Uyw_c + v_Uy_w + tau*Uyw_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
    const double f_Usx = - Uxs_c + v_Ux_s + tau*Uxs_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);
    const double f_Usy = - Uys_c + v_Uy_s + tau*Uys_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c);

    const double delta_Uwx = (f_Uwx*big_Dx-f_Usx*big_B )/big_detx;
    const double delta_Usx = (f_Usx*big_A -f_Uwx*big_Cx)/big_detx;

    const double delta_Uwy = (f_Uwy*big_Dy-f_Usy*big_B )/big_dety;
    const double delta_Usy = (f_Usy*big_A -f_Uwy*big_Cy)/big_dety;

    error = std::sqrt(delta_Uwx*delta_Uwx + delta_Uwy*delta_Uwy + delta_Usx*delta_Usx + delta_Usy*delta_Usy);

    Uxw_c += delta_Uwx;
    Uyw_c += delta_Uwy;
    Uxs_c += delta_Usx;
    Uys_c += delta_Usy;

    //std::cout << count << " " << error << std::endl;

  }


}



void
TG2_scheme::compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant)
{

  // look at tmesh.h 
  const auto & index_quadrant = quadrant->get_global_quad_idx (); 
  const auto & index_quadrant_local = quadrant->get_forest_quad_idx ();

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
      hwdof [ii] = sol [ordhw  (quadrant->gt (ii) )];
      hsdof [ii] = sol [ordhs  (quadrant->gt (ii) )];
      Uxwdof[ii] = sol [ordUxw (quadrant->gt (ii) )];
      Uywdof[ii] = sol [ordUyw (quadrant->gt (ii) )]; 
      Uxsdof[ii] = sol [ordUxs (quadrant->gt (ii) )];
      Uysdof[ii] = sol [ordUys (quadrant->gt (ii) )]; 

      Z_node[ii]  = Z [quadrant->gt (ii)];
      
      isdof_or_hanging[ii] = 1.;
    }
    else
    {
      hwdof [ii] = .5 * (sol [ordhw  (quadrant->gparent (0, ii) )] +
                         sol [ordhw  (quadrant->gparent (1, ii) )]);
      hsdof [ii] = .5 * (sol [ordhs  (quadrant->gparent (0, ii) )] +
                         sol [ordhs  (quadrant->gparent (1, ii) )]);
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
    const auto lambdas = max_eigen (hwdof[ii], hsdof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii]);
    vel_rusanov_cell_x += lambdas[0];
    vel_rusanov_cell_y += lambdas[1];
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.;

  const std::array<double,4> eta_vec = {hdof[0]+Z_node[0], hdof[1]+Z_node[1], hdof[2]+Z_node[2], hdof[3]+Z_node[3]};
 
  grad_cell_hw    = {.5 * ( (hwdof   [3] - hwdof   [2]) + (hwdof   [1] - hwdof   [0]) ), .5 * ( (hwdof   [2] - hwdof   [0]) + (hwdof   [3] - hwdof   [1]) )};
  grad_cell_Uxw   = {.5 * ( (Uxwdof  [3] - Uxwdof  [2]) + (Uxwdof  [1] - Uxwdof  [0]) ), .5 * ( (Uxwdof  [2] - Uxwdof  [0]) + (Uxwdof  [3] - Uxwdof  [1]) )};
  grad_cell_Uyw   = {.5 * ( (Uywdof  [3] - Uywdof  [2]) + (Uywdof  [1] - Uywdof  [0]) ), .5 * ( (Uywdof  [2] - Uywdof  [0]) + (Uywdof  [3] - Uywdof  [1]) )};

  grad_cell_hs    = {.5 * ( (hsdof   [3] - hsdof   [2]) + (hsdof   [1] - hsdof   [0]) ), .5 * ( (hsdof   [2] - hsdof   [0]) + (hsdof   [3] - hsdof   [1]) )};
  grad_cell_Uxs   = {.5 * ( (Uxsdof  [3] - Uxsdof  [2]) + (Uxsdof  [1] - Uxsdof  [0]) ), .5 * ( (Uxsdof  [2] - Uxsdof  [0]) + (Uxsdof  [3] - Uxsdof  [1]) )};
  grad_cell_Uys   = {.5 * ( (Uysdof  [3] - Uysdof  [2]) + (Uysdof  [1] - Uysdof  [0]) ), .5 * ( (Uysdof  [2] - Uysdof  [0]) + (Uysdof  [3] - Uysdof  [1]) )};


  const double & hw_cell    = sol_onehalf[ordhw    (index_quadrant)];
  const double & Uxw_cell   = sol_onehalf[ordUxw   (index_quadrant)];
  const double & Uyw_cell   = sol_onehalf[ordUyw   (index_quadrant)];

  const double & hs_cell    = sol_onehalf[ordhs    (index_quadrant)];
  const double & Uxs_cell   = sol_onehalf[ordUxs   (index_quadrant)];
  const double & Uys_cell   = sol_onehalf[ordUys   (index_quadrant)];

  const double & Z_cell     = Z_onehalf[index_quadrant];

  const auto slope_x_cell = slope_x[index_quadrant_local];
  const auto slope_y_cell = slope_y[index_quadrant_local];


  const auto contr_slope_xs = .5*.5*src_slope_formula (hs_cell, slope_x_cell)*Dx*Dy;
  const auto contr_slope_ys = .5*.5*src_slope_formula (hs_cell, slope_y_cell)*Dx*Dy;

  const auto contr_slope_xw = .5*.5*src_slope_formula (hw_cell, slope_x_cell)*Dx*Dy;
  const auto contr_slope_yw = .5*.5*src_slope_formula (hw_cell, slope_y_cell)*Dx*Dy;


  const auto h_cell  = hw_cell+hs_cell;
  const auto n_cell  = h_cell>epsilon ? hw_cell/h_cell : 0.;
  const auto ns_cell = h_cell>epsilon ? hs_cell/h_cell : 0.;

  std::array<double,4> contr_x_w = {0., 0., 0., 0.}, contr_y_w = {0., 0., 0., 0.}, contr_x_s = {0., 0., 0., 0.}, contr_y_s = {0., 0., 0., 0.};



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

      const auto & index_quadrant_nei = quadrant_nei->get_global_quad_idx ();


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
          double Z_cell_nei, hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei;

          Z_cell_nei     = Z_onehalf  [          index_quadrant_nei ];
          hw_cell_nei    = sol_onehalf[ordhw    (index_quadrant_nei)];
          hs_cell_nei    = sol_onehalf[ordhs    (index_quadrant_nei)];
          Uxw_cell_nei   = sol_onehalf[ordUxw   (index_quadrant_nei)];
          Uyw_cell_nei   = sol_onehalf[ordUyw   (index_quadrant_nei)];
          Uxs_cell_nei   = sol_onehalf[ordUxs   (index_quadrant_nei)];
          Uys_cell_nei   = sol_onehalf[ordUys   (index_quadrant_nei)];

          // .5 salta fuori dall'integrazione per trapezi tra 0 e 1 in coordinata \xi (è il valore in LHS da metter qui sotto!)
          contr_x_w[i_1] += .5*signum(outward_normal_edge[0])*(grav*hw_cell*(Z_cell_nei - Z_cell) +         grav*hw_cell*(hs_cell_nei - hs_cell))*isdof_or_hanging[i_1];
          contr_x_w[i_2] += .5*signum(outward_normal_edge[0])*(grav*hw_cell*(Z_cell_nei - Z_cell) +         grav*hw_cell*(hs_cell_nei - hs_cell))*isdof_or_hanging[i_2];

          contr_y_w[i_1] += .5*signum(outward_normal_edge[1])*(grav*hw_cell*(Z_cell_nei - Z_cell) +         grav*hw_cell*(hs_cell_nei - hs_cell))*isdof_or_hanging[i_1];
          contr_y_w[i_2] += .5*signum(outward_normal_edge[1])*(grav*hw_cell*(Z_cell_nei - Z_cell) +         grav*hw_cell*(hs_cell_nei - hs_cell))*isdof_or_hanging[i_2];


          contr_x_s[i_1] += .5*signum(outward_normal_edge[0])*(grav*hs_cell*(Z_cell_nei - Z_cell) + r_coeff*grav*hs_cell*(hw_cell_nei - hw_cell))*isdof_or_hanging[i_1];
          contr_x_s[i_2] += .5*signum(outward_normal_edge[0])*(grav*hs_cell*(Z_cell_nei - Z_cell) + r_coeff*grav*hs_cell*(hw_cell_nei - hw_cell))*isdof_or_hanging[i_2];

          contr_y_s[i_1] += .5*signum(outward_normal_edge[1])*(grav*hs_cell*(Z_cell_nei - Z_cell) + r_coeff*grav*hs_cell*(hw_cell_nei - hw_cell))*isdof_or_hanging[i_1];
          contr_y_s[i_2] += .5*signum(outward_normal_edge[1])*(grav*hs_cell*(Z_cell_nei - Z_cell) + r_coeff*grav*hs_cell*(hw_cell_nei - hw_cell))*isdof_or_hanging[i_2];


          //break; // this just goes outside the jEdge cycle 
        }
      }

    }

    if (is_boundary_edge) // set boundary conditions
    { 

      auto hw_cell_nei  = hw_cell;
      auto Uxw_cell_nei = Uxw_cell;
      auto Uyw_cell_nei = Uyw_cell;

      auto hs_cell_nei  = hs_cell;
      auto Uxs_cell_nei = Uxs_cell;
      auto Uys_cell_nei = Uys_cell;

      Uxw_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxw_cell + outward_normal_edge[1]*Uyw_cell)*outward_normal_edge[0];
      Uyw_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxw_cell + outward_normal_edge[1]*Uyw_cell)*outward_normal_edge[1];

      Uxs_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxs_cell + outward_normal_edge[1]*Uys_cell)*outward_normal_edge[0];
      Uys_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Uxs_cell + outward_normal_edge[1]*Uys_cell)*outward_normal_edge[1];

      const auto speed     = max_eigen(hw_cell,     hs_cell,     Uxw_cell,     Uyw_cell,     Uxs_cell,     Uys_cell    );
      const auto speed_nei = max_eigen(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei);

      const auto smax = std::max(speed[0]*outward_normal_edge[0]+speed[1]*outward_normal_edge[1], speed_nei[0]*outward_normal_edge[0]+speed_nei[1]*outward_normal_edge[1]); 

      const auto flux_int_hw  = .5*((hw_flux_formula_x (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+hw_flux_formula_x (hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (hw_flux_formula_y (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+hw_flux_formula_y (hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(hw_cell_nei -hw_cell );
      const auto flux_int_Uxw = .5*((Uxw_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uxw_flux_formula_x(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (Uxw_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uxw_flux_formula_y(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uxw_cell_nei-Uxw_cell);
      const auto flux_int_Uyw = .5*((Uyw_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uyw_flux_formula_x(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (Uyw_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uyw_flux_formula_y(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uyw_cell_nei-Uyw_cell);

      const auto flux_int_hs  = .5*((hs_flux_formula_x (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+hs_flux_formula_x (hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (hs_flux_formula_y (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+hs_flux_formula_y (hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(hs_cell_nei -hs_cell );
      const auto flux_int_Uxs = .5*((Uxs_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uxs_flux_formula_x(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (Uxs_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uxs_flux_formula_y(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uxs_cell_nei-Uxs_cell);
      const auto flux_int_Uys = .5*((Uys_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uys_flux_formula_x(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[0] + (Uys_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell)+Uys_flux_formula_y(hw_cell_nei, hs_cell_nei, Uxw_cell_nei, Uyw_cell_nei, Uxs_cell_nei, Uys_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uys_cell_nei-Uys_cell);

      
      // .5 is the base function evaluated in the middle, mid-point intergration
      if (! quadrant->is_hanging (i_1))
      {
        incr[ordhw   (quadrant->gt (i_1))] += -edge_length*flux_int_hw *.5;
        incr[ordUxw  (quadrant->gt (i_1))] += -edge_length*flux_int_Uxw*.5 + contr_slope_xw*outward_normal_edge[0];
        incr[ordUyw  (quadrant->gt (i_1))] += -edge_length*flux_int_Uyw*.5 + contr_slope_yw*outward_normal_edge[1];

        incr[ordhs   (quadrant->gt (i_1))] += -edge_length*flux_int_hs *.5;
        incr[ordUxs  (quadrant->gt (i_1))] += -edge_length*flux_int_Uxs*.5 + contr_slope_xs*outward_normal_edge[0];
        incr[ordUys  (quadrant->gt (i_1))] += -edge_length*flux_int_Uys*.5 + contr_slope_ys*outward_normal_edge[1];
      }
      else
      {
        incr [ordhw  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_hw *.5*.5;
        incr [ordhw  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_hw *.5*.5;
      
        incr [ordUxw (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uxw*.5*.5 + contr_slope_xw*.5*outward_normal_edge[0];
        incr [ordUxw (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uxw*.5*.5 + contr_slope_xw*.5*outward_normal_edge[0];
      
        incr [ordUyw (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uyw*.5*.5 + contr_slope_yw*.5*outward_normal_edge[1];
        incr [ordUyw (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uyw*.5*.5 + contr_slope_yw*.5*outward_normal_edge[1];


        incr [ordhs  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_hs *.5*.5;
        incr [ordhs  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_hs *.5*.5;
      
        incr [ordUxs (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uxs*.5*.5 + contr_slope_xs*.5*outward_normal_edge[0];
        incr [ordUxs (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uxs*.5*.5 + contr_slope_xs*.5*outward_normal_edge[0];
      
        incr [ordUys (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uys*.5*.5 + contr_slope_ys*.5*outward_normal_edge[1];
        incr [ordUys (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uys*.5*.5 + contr_slope_ys*.5*outward_normal_edge[1];
      }

      if (! quadrant->is_hanging (i_2))
      {
        incr[ordhw   (quadrant->gt (i_2))] += -edge_length*flux_int_hw *.5;
        incr[ordUxw  (quadrant->gt (i_2))] += -edge_length*flux_int_Uxw*.5 + contr_slope_xw*outward_normal_edge[0];
        incr[ordUyw  (quadrant->gt (i_2))] += -edge_length*flux_int_Uyw*.5 + contr_slope_yw*outward_normal_edge[1];

        incr[ordhs   (quadrant->gt (i_2))] += -edge_length*flux_int_hs *.5;
        incr[ordUxs  (quadrant->gt (i_2))] += -edge_length*flux_int_Uxs*.5 + contr_slope_xs*outward_normal_edge[0];
        incr[ordUys  (quadrant->gt (i_2))] += -edge_length*flux_int_Uys*.5 + contr_slope_ys*outward_normal_edge[1];
      }
      else
      {
        // il secondo .5 è per hanging nodes
        incr [ordhw  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_hw *.5*.5;
        incr [ordhw  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_hw *.5*.5;
      
        incr [ordUxw (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uxw*.5*.5 + contr_slope_xw*.5*outward_normal_edge[0];
        incr [ordUxw (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uxw*.5*.5 + contr_slope_xw*.5*outward_normal_edge[0];
      
        incr [ordUyw (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uyw*.5*.5 + contr_slope_yw*.5*outward_normal_edge[1];
        incr [ordUyw (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uyw*.5*.5 + contr_slope_yw*.5*outward_normal_edge[1];


        incr [ordhs  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_hs *.5*.5;
        incr [ordhs  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_hs *.5*.5;
      
        incr [ordUxs (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uxs*.5*.5 + contr_slope_xs*.5*outward_normal_edge[0];
        incr [ordUxs (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uxs*.5*.5 + contr_slope_xs*.5*outward_normal_edge[0];
      
        incr [ordUys (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uys*.5*.5 + contr_slope_ys*.5*outward_normal_edge[1];
        incr [ordUys (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uys*.5*.5 + contr_slope_ys*.5*outward_normal_edge[1];
      }

    }

  }
  

  const auto diff_term_hw_x  = h_cell>epsilon ? n_cell*grad_cell_eta[0]*vel_rusanov_cell_y*.5 : grad_cell_hw[0]*vel_rusanov_cell_y*.5;
  const auto diff_term_hw_y  = h_cell>epsilon ? n_cell*grad_cell_eta[1]*vel_rusanov_cell_x*.5 : grad_cell_hw[1]*vel_rusanov_cell_x*.5;

  const auto diff_term_hs_x  = h_cell>epsilon ? ns_cell*grad_cell_eta[1]*vel_rusanov_cell_x*.5 : grad_cell_hs[0]*vel_rusanov_cell_y*.5;
  const auto diff_term_hs_y  = h_cell>epsilon ? ns_cell*grad_cell_eta[1]*vel_rusanov_cell_x*.5 : grad_cell_hs[1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uxw_x = grad_cell_Uxw [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uxw_y = grad_cell_Uxw [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uyw_x = grad_cell_Uyw [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uyw_y = grad_cell_Uyw [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uxs_x = grad_cell_Uxs [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uxs_y = grad_cell_Uxs [1]*vel_rusanov_cell_x*.5;

  const auto diff_term_Uys_x = grad_cell_Uys [0]*vel_rusanov_cell_y*.5;
  const auto diff_term_Uys_y = grad_cell_Uys [1]*vel_rusanov_cell_x*.5;




  const auto F_star_hw_x  = hw_flux_formula_x (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_hw_x;
  const auto F_star_hw_y  = hw_flux_formula_y (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_hw_y;

  const auto F_star_hs_x  = hs_flux_formula_x (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_hs_x;
  const auto F_star_hs_y  = hs_flux_formula_y (hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_hs_y;

  const auto F_star_Uxw_x = Uxw_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uxw_x;
  const auto F_star_Uxw_y = Uxw_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uxw_y;

  const auto F_star_Uyw_x = Uyw_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uyw_x;
  const auto F_star_Uyw_y = Uyw_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uyw_y;

  const auto F_star_Uxs_x = Uxs_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uxs_x;
  const auto F_star_Uxs_y = Uxs_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uxs_y;

  const auto F_star_Uys_x = Uys_flux_formula_x(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uys_x;
  const auto F_star_Uys_y = Uys_flux_formula_y(hw_cell, hs_cell, Uxw_cell, Uyw_cell, Uxs_cell, Uys_cell) - diff_term_Uys_y;



  for (int ii = 0; ii < 4; ++ii){

    const auto hw_    = der_coeffs_x[ii]*F_star_hw_x +der_coeffs_y[ii]*F_star_hw_y;
    const auto Uxw_   = der_coeffs_x[ii]*F_star_Uxw_x+der_coeffs_y[ii]*F_star_Uxw_y - .5*Dy*contr_x_w[ii];
    const auto Uyw_   = der_coeffs_x[ii]*F_star_Uyw_x+der_coeffs_y[ii]*F_star_Uyw_y - .5*Dx*contr_y_w[ii];

    const auto hs_    = der_coeffs_x[ii]*F_star_hs_x +der_coeffs_y[ii]*F_star_hs_y;
    const auto Uxs_   = der_coeffs_x[ii]*F_star_Uxs_x+der_coeffs_y[ii]*F_star_Uxs_y - .5*Dy*contr_x_s[ii];
    const auto Uys_   = der_coeffs_x[ii]*F_star_Uys_x+der_coeffs_y[ii]*F_star_Uys_y - .5*Dx*contr_y_s[ii];


    const auto hw_al  = der_coeffs_x[ii]*diff_term_hw_x  + der_coeffs_y[ii]*diff_term_hw_y; 
    const auto Uxw_al = der_coeffs_x[ii]*diff_term_Uxw_x + der_coeffs_y[ii]*diff_term_Uxw_y;
    const auto Uyw_al = der_coeffs_x[ii]*diff_term_Uyw_x + der_coeffs_y[ii]*diff_term_Uyw_y;

    const auto hs_al  = der_coeffs_x[ii]*diff_term_hs_x  + der_coeffs_y[ii]*diff_term_hs_y; 
    const auto Uxs_al = der_coeffs_x[ii]*diff_term_Uxs_x + der_coeffs_y[ii]*diff_term_Uxs_y;
    const auto Uys_al = der_coeffs_x[ii]*diff_term_Uys_x + der_coeffs_y[ii]*diff_term_Uys_y;

    incr_anti_diff[ordhw (index_quadrant_local)][ii] = hw_al;
    incr_anti_diff[ordUxw(index_quadrant_local)][ii] = Uxw_al;
    incr_anti_diff[ordUyw(index_quadrant_local)][ii] = Uyw_al;

    incr_anti_diff[ordhs (index_quadrant_local)][ii] = hs_al;
    incr_anti_diff[ordUxs(index_quadrant_local)][ii] = Uxs_al;
    incr_anti_diff[ordUys(index_quadrant_local)][ii] = Uys_al;


    if (! quadrant->is_hanging (ii)){

      incr [ordhw  (quadrant->gt (ii))] += hw_;
      incr [ordUxw (quadrant->gt (ii))] += Uxw_;
      incr [ordUyw (quadrant->gt (ii))] += Uyw_;

      P_plus [ordhw  (quadrant->gt (ii))] += std::max(0., hw_al );
      P_plus [ordUxw (quadrant->gt (ii))] += std::max(0., Uxw_al);
      P_plus [ordUyw (quadrant->gt (ii))] += std::max(0., Uyw_al);

      P_minus [ordhw  (quadrant->gt (ii))] += std::min(0., hw_al );
      P_minus [ordUxw (quadrant->gt (ii))] += std::min(0., Uxw_al);
      P_minus [ordUyw (quadrant->gt (ii))] += std::min(0., Uyw_al);



      incr [ordhs  (quadrant->gt (ii))] += hs_;
      incr [ordUxs (quadrant->gt (ii))] += Uxs_;
      incr [ordUys (quadrant->gt (ii))] += Uys_;

      P_plus [ordhs  (quadrant->gt (ii))] += std::max(0., hs_al );
      P_plus [ordUxs (quadrant->gt (ii))] += std::max(0., Uxs_al);
      P_plus [ordUys (quadrant->gt (ii))] += std::max(0., Uys_al);

      P_minus [ordhs  (quadrant->gt (ii))] += std::min(0., hs_al );
      P_minus [ordUxs (quadrant->gt (ii))] += std::min(0., Uxs_al);
      P_minus [ordUys (quadrant->gt (ii))] += std::min(0., Uys_al);
      
      
    } else {

      // w part
      incr [ordhw  (quadrant->gparent(0,ii))] += hw_;
      incr [ordhw  (quadrant->gparent(1,ii))] += hw_;
      
      incr [ordUxw (quadrant->gparent(0,ii))] += Uxw_;
      incr [ordUxw (quadrant->gparent(1,ii))] += Uxw_;
      
      incr [ordUyw (quadrant->gparent(0,ii))] += Uyw_;
      incr [ordUyw (quadrant->gparent(1,ii))] += Uyw_;



      P_plus [ordhw  (quadrant->gparent(0,ii))] += std::max(0., hw_al);
      P_plus [ordhw  (quadrant->gparent(1,ii))] += std::max(0., hw_al);
      
      P_plus [ordUxw (quadrant->gparent(0,ii))] += std::max(0., Uxw_al);
      P_plus [ordUxw (quadrant->gparent(1,ii))] += std::max(0., Uxw_al);
      
      P_plus [ordUyw (quadrant->gparent(0,ii))] += std::max(0., Uyw_al);
      P_plus [ordUyw (quadrant->gparent(1,ii))] += std::max(0., Uyw_al);



      P_minus [ordhw  (quadrant->gparent(0,ii))] += std::min(0., hw_al);
      P_minus [ordhw  (quadrant->gparent(1,ii))] += std::min(0., hw_al);
      
      P_minus [ordUxw (quadrant->gparent(0,ii))] += std::min(0., Uxw_al);
      P_minus [ordUxw (quadrant->gparent(1,ii))] += std::min(0., Uxw_al);
      
      P_minus [ordUyw (quadrant->gparent(0,ii))] += std::min(0., Uyw_al);
      P_minus [ordUyw (quadrant->gparent(1,ii))] += std::min(0., Uyw_al);



      // s part
      incr [ordhs  (quadrant->gparent(0,ii))] += hs_;
      incr [ordhs  (quadrant->gparent(1,ii))] += hs_;
      
      incr [ordUxs (quadrant->gparent(0,ii))] += Uxs_;
      incr [ordUxs (quadrant->gparent(1,ii))] += Uxs_;
      
      incr [ordUys (quadrant->gparent(0,ii))] += Uys_;
      incr [ordUys (quadrant->gparent(1,ii))] += Uys_;



      P_plus [ordhs  (quadrant->gparent(0,ii))] += std::max(0., hs_al);
      P_plus [ordhs  (quadrant->gparent(1,ii))] += std::max(0., hs_al);
      
      P_plus [ordUxs (quadrant->gparent(0,ii))] += std::max(0., Uxs_al);
      P_plus [ordUxs (quadrant->gparent(1,ii))] += std::max(0., Uxs_al);
      
      P_plus [ordUys (quadrant->gparent(0,ii))] += std::max(0., Uys_al);
      P_plus [ordUys (quadrant->gparent(1,ii))] += std::max(0., Uys_al);



      P_minus [ordhs  (quadrant->gparent(0,ii))] += std::min(0., hs_al);
      P_minus [ordhs  (quadrant->gparent(1,ii))] += std::min(0., hs_al);
      
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

    double hwdof_c, Uxwdof_c, Uywdof_c, P_plus_hw_c, P_minus_hw_c, P_plus_Uxw_c, P_minus_Uxw_c, P_plus_Uyw_c, P_minus_Uyw_c;
    double hsdof_c, Uxsdof_c, Uysdof_c, P_plus_hs_c, P_minus_hs_c, P_plus_Uxs_c, P_minus_Uxs_c, P_plus_Uys_c, P_minus_Uys_c;

    if (! quadrant->is_hanging (ii)){
      hwdof_c      = sol [ordhw    (quadrant->gt (ii))];
      Uxwdof_c     = sol [ordUxw   (quadrant->gt (ii))];
      Uywdof_c     = sol [ordUyw   (quadrant->gt (ii))];

      hsdof_c      = sol [ordhs    (quadrant->gt (ii))];
      Uxsdof_c     = sol [ordUxs   (quadrant->gt (ii))];
      Uysdof_c     = sol [ordUys   (quadrant->gt (ii))];

      Z_node[ii]    = Z [quadrant->gt (ii)];

      P_plus_hw_c   = P_plus [ordhw   (quadrant->gt (ii))];
      P_minus_hw_c  = P_minus[ordhw   (quadrant->gt (ii))];

      P_plus_Uxw_c   = P_plus [ordUxw   (quadrant->gt (ii))];
      P_minus_Uxw_c  = P_minus[ordUxw   (quadrant->gt (ii))];

      P_plus_Uyw_c   = P_plus [ordUyw   (quadrant->gt (ii))];
      P_minus_Uyw_c  = P_minus[ordUyw   (quadrant->gt (ii))];


      P_plus_hs_c   = P_plus [ordhs   (quadrant->gt (ii))];
      P_minus_hs_c  = P_minus[ordhs   (quadrant->gt (ii))];

      P_plus_Uxs_c   = P_plus [ordUxs   (quadrant->gt (ii))];
      P_minus_Uxs_c  = P_minus[ordUxs   (quadrant->gt (ii))];

      P_plus_Uys_c   = P_plus [ordUys   (quadrant->gt (ii))];
      P_minus_Uys_c  = P_minus[ordUys   (quadrant->gt (ii))];
      

    } else {
      hwdof_c   = .5 * (sol [ordhw  (quadrant->gparent(0,ii))] +
                        sol [ordhw  (quadrant->gparent(1,ii))]);
      Uxwdof_c  = .5 * (sol [ordUxw (quadrant->gparent(0,ii))] +
                        sol [ordUxw (quadrant->gparent(1,ii))]);
      Uywdof_c  = .5 * (sol [ordUyw (quadrant->gparent(0,ii))] +
                        sol [ordUyw (quadrant->gparent(1,ii))]);


      hsdof_c   = .5 * (sol [ordhs  (quadrant->gparent(0,ii))] +
                        sol [ordhs  (quadrant->gparent(1,ii))]);
      Uxsdof_c  = .5 * (sol [ordUxs (quadrant->gparent(0,ii))] +
                        sol [ordUxs (quadrant->gparent(1,ii))]);
      Uysdof_c  = .5 * (sol [ordUys (quadrant->gparent(0,ii))] +
                        sol [ordUys (quadrant->gparent(1,ii))]);

      Z_node[ii] = .5 * (Z [quadrant->gparent(0,ii)] +
                         Z [quadrant->gparent(1,ii)]);

      P_plus_hw_c   = .5 * (P_plus [ordhw (quadrant->gparent(0,ii))] +
                            P_plus [ordhw (quadrant->gparent(1,ii))]);
      P_minus_hw_c  = .5 * (P_minus [ordhw (quadrant->gparent(0,ii))] +
                            P_minus [ordhw (quadrant->gparent(1,ii))]);

      P_plus_Uxw_c   = .5 * (P_plus [ordUxw (quadrant->gparent(0,ii))] +
                             P_plus [ordUxw (quadrant->gparent(1,ii))]);
      P_minus_Uxw_c  = .5 * (P_minus [ordUxw (quadrant->gparent(0,ii))] +
                             P_minus [ordUxw (quadrant->gparent(1,ii))]);

      P_plus_Uyw_c   = .5 * (P_plus [ordUyw (quadrant->gparent(0,ii))] +
                             P_plus [ordUyw (quadrant->gparent(1,ii))]);
      P_minus_Uyw_c  = .5 * (P_minus [ordUyw (quadrant->gparent(0,ii))] +
                             P_minus [ordUyw (quadrant->gparent(1,ii))]);



      P_plus_hs_c   = .5 * (P_plus [ordhs (quadrant->gparent(0,ii))] +
                            P_plus [ordhs (quadrant->gparent(1,ii))]);
      P_minus_hs_c  = .5 * (P_minus [ordhs (quadrant->gparent(0,ii))] +
                            P_minus [ordhs (quadrant->gparent(1,ii))]);

      P_plus_Uxs_c   = .5 * (P_plus [ordUxs (quadrant->gparent(0,ii))] +
                             P_plus [ordUxs (quadrant->gparent(1,ii))]);
      P_minus_Uxs_c  = .5 * (P_minus [ordUxs (quadrant->gparent(0,ii))] +
                             P_minus [ordUxs (quadrant->gparent(1,ii))]);

      P_plus_Uys_c   = .5 * (P_plus [ordUys (quadrant->gparent(0,ii))] +
                             P_plus [ordUys (quadrant->gparent(1,ii))]);
      P_minus_Uys_c  = .5 * (P_minus [ordUys (quadrant->gparent(0,ii))] +
                             P_minus [ordUys (quadrant->gparent(1,ii))]);
      
    }


    double hdof_c  = hwdof_c+hsdof_c;
    double ndof_c  = hdof_c>epsilon ? hwdof_c/hdof_c : 0.;
    double nsdof_c = hdof_c>epsilon ? hsdof_c/hdof_c : 0.;

    hwdof      [ii] = hdof_c>epsilon ? hwdof_c+Z_node[ii]*ndof_c  : hwdof_c;
    hsdof      [ii] = hdof_c>epsilon ? hsdof_c+Z_node[ii]*nsdof_c : hsdof_c;
    Uxwdof     [ii] = Uxwdof_c;
    Uywdof     [ii] = Uywdof_c;
    Uxsdof     [ii] = Uxsdof_c;
    Uysdof     [ii] = Uysdof_c;

    P_plus_hw_dof [ii] = P_plus_hw_c;
    P_minus_hw_dof[ii] = P_minus_hw_c;

    P_plus_hs_dof [ii] = P_plus_hs_c;
    P_minus_hs_dof[ii] = P_minus_hs_c;

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
  const auto hw_min_cell  = *std::min_element(hwdof.begin() ,  hwdof.end() );
  const auto hw_max_cell  = *std::max_element(hwdof.begin() ,  hwdof.end() );

  const auto hs_min_cell  = *std::min_element(hsdof.begin() ,  hsdof.end() );
  const auto hs_max_cell  = *std::max_element(hsdof.begin() ,  hsdof.end() );
 
  const auto Uxw_min_cell = *std::min_element(Uxwdof.begin(),  Uxwdof.end());
  const auto Uxw_max_cell = *std::max_element(Uxwdof.begin(),  Uxwdof.end());

  const auto Uyw_min_cell = *std::min_element(Uywdof.begin(),  Uywdof.end());
  const auto Uyw_max_cell = *std::max_element(Uywdof.begin(),  Uywdof.end());

  const auto Uxs_min_cell = *std::min_element(Uxsdof.begin(),  Uxsdof.end());
  const auto Uxs_max_cell = *std::max_element(Uxsdof.begin(),  Uxsdof.end());

  const auto Uys_min_cell = *std::min_element(Uysdof.begin(),  Uysdof.end());
  const auto Uys_max_cell = *std::max_element(Uysdof.begin(),  Uysdof.end());   

  bool is_node_in_element = false;

  std::array<double,4> hw_min  = {hw_min_cell, hw_min_cell, hw_min_cell, hw_min_cell }, hw_max  = {hw_max_cell, hw_max_cell, hw_max_cell, hw_max_cell },
                       hs_min  = {hs_min_cell, hs_min_cell, hs_min_cell, hs_min_cell }, hs_max  = {hs_max_cell, hs_max_cell, hs_max_cell, hs_max_cell },
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

          double n_current_cell, ns_current_cell, hw_current_cell, hs_current_cell, Uxw_current_cell, Uyw_current_cell, Uxs_current_cell, Uys_current_cell;

          if (! quadrant_nei->is_hanging (jj)){
            n_current_cell  = (sol [ordhw  (quadrant_nei->gt (jj))]+sol [ordhs  (quadrant_nei->gt (jj))])>epsilon ? sol [ordhw  (quadrant_nei->gt (jj))]/(sol [ordhw  (quadrant_nei->gt (jj))]+sol [ordhs  (quadrant_nei->gt (jj))]) : 0.;
            ns_current_cell = (sol [ordhw  (quadrant_nei->gt (jj))]+sol [ordhs  (quadrant_nei->gt (jj))])>epsilon ? sol [ordhs  (quadrant_nei->gt (jj))]/(sol [ordhw  (quadrant_nei->gt (jj))]+sol [ordhs  (quadrant_nei->gt (jj))]) : 0.;

            hw_current_cell  = sol [ordhw  (quadrant_nei->gt (jj))] + Z [quadrant_nei->gt (jj)]*n_current_cell;
            hs_current_cell  = sol [ordhs  (quadrant_nei->gt (jj))] + Z [quadrant_nei->gt (jj)]*ns_current_cell;
            Uxw_current_cell = sol [ordUxw (quadrant_nei->gt (jj))];
            Uyw_current_cell = sol [ordUyw (quadrant_nei->gt (jj))];
            Uxs_current_cell = sol [ordUxs (quadrant_nei->gt (jj))];
            Uys_current_cell = sol [ordUys (quadrant_nei->gt (jj))];
          } else {
            
            n_current_cell = (.5 * (sol [ordhw  (quadrant_nei->gparent(0,jj))] +
                                    sol [ordhw  (quadrant_nei->gparent(1,jj))])+
                              .5 * (sol [ordhs  (quadrant_nei->gparent(0,jj))] +
                                    sol [ordhs  (quadrant_nei->gparent(1,jj))]))>epsilon ? (.5 * (sol [ordhw  (quadrant_nei->gparent(0,jj))] + sol [ordhw  (quadrant_nei->gparent(1,jj))]))/(.5 * (sol [ordhw  (quadrant_nei->gparent(0,jj))]+sol [ordhw  (quadrant_nei->gparent(1,jj))])+.5 * (sol [ordhs  (quadrant_nei->gparent(0,jj))] + sol [ordhs  (quadrant_nei->gparent(1,jj))])) : 0.;

            ns_current_cell = (.5 * (sol [ordhw  (quadrant_nei->gparent(0,jj))] +
                                    sol [ordhw  (quadrant_nei->gparent(1,jj))])+
                               .5 * (sol [ordhs  (quadrant_nei->gparent(0,jj))] +
                                    sol [ordhs  (quadrant_nei->gparent(1,jj))]))>epsilon ? (.5 * (sol [ordhs  (quadrant_nei->gparent(0,jj))] + sol [ordhs  (quadrant_nei->gparent(1,jj))]))/(.5 * (sol [ordhw  (quadrant_nei->gparent(0,jj))]+sol [ordhw  (quadrant_nei->gparent(1,jj))])+.5 * (sol [ordhs  (quadrant_nei->gparent(0,jj))] + sol [ordhs  (quadrant_nei->gparent(1,jj))])) : 0.;


            hw_current_cell  = .5 * (sol [ordhw  (quadrant_nei->gparent(0,jj))] +
                                     sol [ordhw  (quadrant_nei->gparent(1,jj))]) + 
                               .5 * (Z [quadrant_nei->gparent(0,jj)] +
                                     Z [quadrant_nei->gparent(1,jj)])*n_current_cell;
            hs_current_cell  = .5 * (sol [ordhs  (quadrant_nei->gparent(0,jj))] +
                                     sol [ordhs  (quadrant_nei->gparent(1,jj))]) +
                               .5 * (Z [quadrant_nei->gparent(0,jj)] +
                                     Z [quadrant_nei->gparent(1,jj)])*ns_current_cell;
            Uxw_current_cell = .5 * (sol [ordUxw (quadrant_nei->gparent(0,jj))] +
                                     sol [ordUxw (quadrant_nei->gparent(1,jj))]);
            Uyw_current_cell = .5 * (sol [ordUyw (quadrant_nei->gparent(0,jj))] +
                                     sol [ordUyw (quadrant_nei->gparent(1,jj))]);
            Uxs_current_cell = .5 * (sol [ordUxs (quadrant_nei->gparent(0,jj))] +
                                     sol [ordUxs (quadrant_nei->gparent(1,jj))]);
            Uys_current_cell = .5 * (sol [ordUys (quadrant_nei->gparent(0,jj))] +
                                     sol [ordUys (quadrant_nei->gparent(1,jj))]);
          }

          hw_min[ii]  = std::min(hw_min[ii],  hw_current_cell );
          hw_max[ii]  = std::max(hw_max[ii],  hw_current_cell );

          hs_min[ii]  = std::min(hs_min[ii],  hs_current_cell );
          hs_max[ii]  = std::max(hs_max[ii],  hs_current_cell );

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
  double phi_cell_hw = 1., phi_cell_hs = 1., phi_cell_Uxw = 1., phi_cell_Uyw = 1., phi_cell_Uxs = 1., phi_cell_Uys = 1.;
  for (int ii = 0; ii < 4; ++ii){

    const auto & flux_on_the_node_hw  = incr_anti_diff[ordhw (index_quadrant)][ii];
    const auto & flux_on_the_node_hs  = incr_anti_diff[ordhs (index_quadrant)][ii];
    const auto & flux_on_the_node_Uxw = incr_anti_diff[ordUxw(index_quadrant)][ii];
    const auto & flux_on_the_node_Uyw = incr_anti_diff[ordUyw(index_quadrant)][ii];
    const auto & flux_on_the_node_Uxs = incr_anti_diff[ordUxs(index_quadrant)][ii];
    const auto & flux_on_the_node_Uys = incr_anti_diff[ordUys(index_quadrant)][ii];


    const auto speed = max_eigen(hwdof[ii], hsdof[ii], Uxwdof[ii], Uywdof[ii], Uxsdof[ii], Uysdof[ii]);
    const auto vel_rusanov_cell_x = speed[0]; //hpoint>epsilon ? (std::abs(Uxdof[ii]/hpoint)+celerity) : 0.;
    const auto vel_rusanov_cell_y = speed[1]; //hpoint>epsilon ? (std::abs(Uydof[ii]/hpoint)+celerity) : 0.;

    const auto vel_square_rusanov_cell = vel_rusanov_cell_x * vel_rusanov_cell_y;

    flux_limiter(hw_min [ii], hw_max [ii], hwdof  [ii], P_plus_hw_dof [ii], P_minus_hw_dof  [ii], flux_on_the_node_hw,  vel_square_rusanov_cell, phi_cell_hw );
    flux_limiter(hs_min [ii], hs_max [ii], hsdof  [ii], P_plus_hs_dof [ii], P_minus_hs_dof  [ii], flux_on_the_node_hs,  vel_square_rusanov_cell, phi_cell_hs );
    flux_limiter(Uxw_min[ii], Uxw_max[ii], Uxwdof [ii], P_plus_Uxw_dof[ii], P_minus_Uxw_dof [ii], flux_on_the_node_Uxw, vel_square_rusanov_cell, phi_cell_Uxw);
    flux_limiter(Uyw_min[ii], Uyw_max[ii], Uywdof [ii], P_plus_Uyw_dof[ii], P_minus_Uyw_dof [ii], flux_on_the_node_Uyw, vel_square_rusanov_cell, phi_cell_Uyw);
    flux_limiter(Uxs_min[ii], Uxs_max[ii], Uxsdof [ii], P_plus_Uxs_dof[ii], P_minus_Uxs_dof [ii], flux_on_the_node_Uxs, vel_square_rusanov_cell, phi_cell_Uxs);
    flux_limiter(Uys_min[ii], Uys_max[ii], Uysdof [ii], P_plus_Uys_dof[ii], P_minus_Uys_dof [ii], flux_on_the_node_Uys, vel_square_rusanov_cell, phi_cell_Uys);
  }

  //std::cout << phi_cell_h << " " << phi_cell_Ux << " " << phi_cell_Uy << std::endl;

  //phi_cell_hw = 0., phi_cell_hs = 0., phi_cell_Uxw = 0., phi_cell_Uyw = 0., phi_cell_Uxs = 0., phi_cell_Uys = 0.; 
  //std::cout << "we put limiter equal to one" << std::endl;



  for (int ii = 0; ii < 4; ++ii){

    //std::cout << phi_cell_h << " " << phi_cell_Ux << " " << phi_cell_Uy << std::endl;

    const auto flux_on_the_node_hw  = incr_anti_diff[ordhw (index_quadrant)][ii]*phi_cell_hw;
    const auto flux_on_the_node_hs  = incr_anti_diff[ordhs (index_quadrant)][ii]*phi_cell_hs;
    const auto flux_on_the_node_Uxw = incr_anti_diff[ordUxw(index_quadrant)][ii]*phi_cell_Uxw;
    const auto flux_on_the_node_Uyw = incr_anti_diff[ordUyw(index_quadrant)][ii]*phi_cell_Uyw;
    const auto flux_on_the_node_Uxs = incr_anti_diff[ordUxs(index_quadrant)][ii]*phi_cell_Uxs;
    const auto flux_on_the_node_Uys = incr_anti_diff[ordUys(index_quadrant)][ii]*phi_cell_Uys;

    if (! quadrant->is_hanging (ii)){

      incr [ordhw  (quadrant->gt (ii))] += flux_on_the_node_hw;
      incr [ordhs  (quadrant->gt (ii))] += flux_on_the_node_hs;
      incr [ordUxw (quadrant->gt (ii))] += flux_on_the_node_Uxw;
      incr [ordUyw (quadrant->gt (ii))] += flux_on_the_node_Uyw;
      incr [ordUxs (quadrant->gt (ii))] += flux_on_the_node_Uxs;
      incr [ordUys (quadrant->gt (ii))] += flux_on_the_node_Uys;

    } else {

      incr [ordhw  (quadrant->gparent(0,ii))] += flux_on_the_node_hw;
      incr [ordhw  (quadrant->gparent(1,ii))] += flux_on_the_node_hw;

      incr [ordhs  (quadrant->gparent(0,ii))] += flux_on_the_node_hs;
      incr [ordhs  (quadrant->gparent(1,ii))] += flux_on_the_node_hs;
      
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
TG2_scheme::loop_step (const int& kk, const bool& isInitial) 
{

  const auto & h_w  = soldd_rkc.get_owned_data ()[kk  ];
  const auto & h_s  = soldd_rkc.get_owned_data ()[kk+1];
  const auto & Ux_w = soldd_rkc.get_owned_data ()[kk+2];
  const auto & Uy_w = soldd_rkc.get_owned_data ()[kk+3]; 
  const auto & Ux_s = soldd_rkc.get_owned_data ()[kk+4];
  const auto & Uy_s = soldd_rkc.get_owned_data ()[kk+5]; 


  const auto hw_ = hw_src_formula(h_w, h_s, Ux_w, Uy_w, Ux_s, Uy_s);
  const auto hs_ = hs_src_formula(h_w, h_s, Ux_w, Uy_w, Ux_s, Uy_s);

  
  if (isInitial)
  {
    incr_initial_source.get_owned_data ()[kk  ] = hw_;
    incr_initial_source.get_owned_data ()[kk+1] = hs_;  
  }
  else
  {
    incr_source.get_owned_data ()[kk  ] = hw_;
    incr_source.get_owned_data ()[kk+1] = hs_;
  }

}

void
TG2_scheme::rkc(const int& j, const int& s, const int& kk)
{
  // kk is the current owned node
  double v_hw, v_hs;

  if (j == 1)
  {
    v_hw = sol_ini_rkc.get_owned_data ()[kk  ];
    v_hs = sol_ini_rkc.get_owned_data ()[kk+1];
  }
  else
  {
    v_hw = (1. - mu_vect[j] - v_vect[j])*sol_ini_rkc.get_owned_data ()[kk  ] + mu_vect[j]*sold_rkc.get_owned_data ()[kk  ] + v_vect[j]*soldd_rkc.get_owned_data ()[kk  ] +  
    (gamma_tilde_vect[j]*mu_tilde_vect[1]*mu_vect[j]/mu_tilde_vect[j] - (1. - mu_vect[j] - v_vect[j])*mu_tilde_vect[1])*dt*incr_initial_source.get_owned_data ()[kk  ] - v_vect[j]*mu_tilde_vect[1]*dt*incr_source.get_owned_data ()[kk  ];

    v_hs = (1. - mu_vect[j] - v_vect[j])*sol_ini_rkc.get_owned_data ()[kk+1] + mu_vect[j]*sold_rkc.get_owned_data ()[kk+1] + v_vect[j]*soldd_rkc.get_owned_data ()[kk+1] +  
    (gamma_tilde_vect[j]*mu_tilde_vect[1]*mu_vect[j]/mu_tilde_vect[j] - (1. - mu_vect[j] - v_vect[j])*mu_tilde_vect[1])*dt*incr_initial_source.get_owned_data ()[kk+1] - v_vect[j]*mu_tilde_vect[1]*dt*incr_source.get_owned_data ()[kk+1];
  }


  auto & hw_c  = sol.get_owned_data ()[kk  ];
  auto & hs_c  = sol.get_owned_data ()[kk+1];
  auto & Uxw_c = sol.get_owned_data ()[kk+2];
  auto & Uyw_c = sol.get_owned_data ()[kk+3];
  auto & Uxs_c = sol.get_owned_data ()[kk+4];
  auto & Uys_c = sol.get_owned_data ()[kk+5];

  hw_c = v_hw;
  hs_c = v_hs;


  const auto U_tot_x = Uxs_c + Uxw_c;
  const auto U_tot_y = Uys_c + Uyw_c;

  const auto abs_mass_flux = std::sqrt(U_tot_x*U_tot_x + U_tot_y*U_tot_y);

  const auto erosion_contribution = erosion_coefficient*abs_mass_flux;


  const double tolerance = 1.e-4;
  const int Nmax = 1e3;
  int count = -1;
  double error = tolerance + 1;
  while (count++<Nmax && error>tolerance)
  {
    const auto h_c = hw_c+hs_c;
    const auto nw_c = h_c>epsilon ? hw_c/h_c : 0.;
    const auto ns_c = h_c>epsilon ? hs_c/h_c : 0.;
 

    // Newton Jacobian matrix, 
    // diagonalizzazione dei termini extra-diag delle matrici g e f
    const double big_A = 1.-(h_c>epsilon ? mu_tilde_vect[1]*dt*ns_c/h_c*erosion_contribution : 0.);
    const double big_B =     h_c>epsilon ? mu_tilde_vect[1]*dt*nw_c/h_c*erosion_contribution : 0.;
    const double big_C =     h_c>epsilon ? mu_tilde_vect[1]*dt*ns_c/h_c*erosion_contribution : 0.;
    const double big_D = 1.-(h_c>epsilon ? mu_tilde_vect[1]*dt*nw_c/h_c*erosion_contribution : 0.);

    // rhs terms of the Newton matrix,
    const auto f_hw = v_hw + mu_tilde_vect[1]*dt*hw_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - hw_c;
    const auto f_hs = v_hs + mu_tilde_vect[1]*dt*hs_src_formula(hw_c, hs_c, Uxw_c, Uyw_c, Uxs_c, Uys_c) - hs_c;

    const double big_det = big_A*big_D-big_B*big_C;


    const auto delta_hw = (f_hw*big_D-f_hs*big_B)/big_det;
    const auto delta_hs = (f_hs*big_A-f_hw*big_C)/big_det;
    

    error = std::sqrt(delta_hw*delta_hw + delta_hs*delta_hs);

    hw_c += delta_hw;
    hs_c += delta_hs;
  }


}


void
TG2_scheme::stabilization_term (const int& kk)
{
  // apply the corrector step now, it's like an interphase drag
  const auto & hw_c  = sol.get_owned_data ()[kk  ];
  const auto & hs_c  = sol.get_owned_data ()[kk+1];
  const auto & Uxw_c = sol.get_owned_data ()[kk+2];
  const auto & Uyw_c = sol.get_owned_data ()[kk+3];
  const auto & Uxs_c = sol.get_owned_data ()[kk+4];
  const auto & Uys_c = sol.get_owned_data ()[kk+5];

  const double h_c = hw_c+hs_c;

  const double vel_w_x = hw_c>epsilon ? Uxw_c/hw_c : 0.;
  const double vel_s_x = hs_c>epsilon ? Uxs_c/hs_c : 0.;
  const double vel_w_y = hw_c>epsilon ? Uyw_c/hw_c : 0.;
  const double vel_s_y = hs_c>epsilon ? Uys_c/hs_c : 0.;

  const auto n  = h_c>epsilon ? hw_c/h_c : 0.;
  const auto ns = h_c>epsilon ? hs_c/h_c : 0.;

  const double kinematic_speed_wave = std::sqrt(grav*h_c);
  const double beta_coeff  = std::sqrt(.5*n*(1.-r_coeff));
  const double beta_coeff_ = std::sqrt(.5*  (1.-r_coeff));

  const auto density = ns*density_s + n*density_w;

  const auto abs_delta_vel_x = std::abs(vel_w_x-vel_s_x);
  const auto abs_delta_vel_y = std::abs(vel_w_y-vel_s_y);

  const auto hyp_diff_x  = abs_delta_vel_x - 2.*kinematic_speed_wave*beta_coeff;
  const auto hyp_diff_y  = abs_delta_vel_y - 2.*kinematic_speed_wave*beta_coeff;
  const auto hyp_diff_x_ = abs_delta_vel_x - 2.*kinematic_speed_wave;
  const auto hyp_diff_y_ = abs_delta_vel_y - 2.*kinematic_speed_wave;

  const double cx_sgn = std::max( h_c>epsilon && hyp_diff_x_<0 ? 100*hyp_diff_x/density/dt*std::sqrt(n)*ns*density_s*density_w/(2.*kinematic_speed_wave*beta_coeff_) : 0., 0.);
  const double cy_sgn = std::max( h_c>epsilon && hyp_diff_y_<0 ? 100*hyp_diff_y/density/dt*std::sqrt(n)*ns*density_s*density_w/(2.*kinematic_speed_wave*beta_coeff_) : 0., 0.);

  //const double cx_sgn = std::max( h_c>epsilon && hyp_diff_x>0 && hyp_diff_x_<0 ? 1.e7 : 0., 0.);
  //const double cy_sgn = std::max( h_c>epsilon && hyp_diff_y>0 && hyp_diff_y_<0 ? 1.e7 : 0., 0.);


  const double big_Ax = 1. + ((hw_c>epsilon) ? dt*cx_sgn/hw_c*h_c/density_w : 0.);
  const double big_Bx = hs_c>epsilon ? -dt*cx_sgn/hs_c*h_c/density_w : 0.;
  const double big_Cx = hw_c>epsilon ? -dt*cx_sgn/hw_c*h_c/density_s : 0.;
  const double big_Dx = 1. + ((hs_c>epsilon) ? dt*cx_sgn/hs_c*h_c/density_s : 0.);

  const double rhs_1x = Uxw_c;
  const double rhs_2x = Uxs_c;

  const double big_detx = big_Ax*big_Dx-big_Bx*big_Cx; 

  const double big_Ay = 1. + ((hw_c>epsilon) ? dt*cy_sgn/hw_c*h_c/density_w : 0.);
  const double big_By = hs_c>epsilon ? -dt*cy_sgn/hs_c*h_c/density_w : 0.;
  const double big_Cy = hw_c>epsilon ? -dt*cy_sgn/hw_c*h_c/density_s : 0.;
  const double big_Dy = 1. + ((hs_c>epsilon) ? dt*cy_sgn/hs_c*h_c/density_s : 0.);

  const double rhs_1y = Uyw_c;
  const double rhs_2y = Uys_c;

  const double big_dety = big_Ay*big_Dy-big_By*big_Cy;


  sol.get_owned_data ()[kk+2] = (rhs_1x*big_Dx-rhs_2x*big_Bx)/big_detx;
  sol.get_owned_data ()[kk+3] = (rhs_1y*big_Dy-rhs_2y*big_By)/big_dety;
  sol.get_owned_data ()[kk+4] = (rhs_2x*big_Ax-rhs_1x*big_Cx)/big_detx;
  sol.get_owned_data ()[kk+5] = (rhs_2y*big_Ay-rhs_1y*big_Cy)/big_dety;

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
TG2_scheme::hw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys) 
{ 
  // A flux-limiting wetting–drying method for finite-element shallow-water models, with application to the Scheldt Estuary
  return (hw>epsilon ? Uxw : 0.); 
}

double
TG2_scheme::hw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  // A flux-limiting wetting–drying method for finite-element shallow-water models, with application to the Scheldt Estuary
  return (hw>epsilon ? Uyw : 0.); 
}

double
TG2_scheme::hs_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  return (hs>epsilon ? Uxs : 0.); 
}

double
TG2_scheme::hs_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  return (hs>epsilon ? Uys : 0.); 
}

double
TG2_scheme::Uxw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto vel_x = hw>epsilon ? Uxw/hw : 0.;
  return (Uxw*vel_x + grav*hw*hw/2.);
}
 
double
TG2_scheme::Uxw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (hw>epsilon ? Uyw*Uxw/hw : 0.); }

double
TG2_scheme::Uyw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (hw>epsilon ? Uyw*Uxw/hw : 0.); }

double
TG2_scheme::Uyw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto vel_y = hw>epsilon ? Uyw/hw : 0.;
  return (Uyw*vel_y + grav*hw*hw/2.); 
}


double
TG2_scheme::Uxs_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto vel_x = hs>epsilon ? Uxs/hs : 0.;
  return (Uxs*vel_x + grav*hs*hs/2. + grav*(1-r_coeff)*hs*hw/2.);  
}
 
double
TG2_scheme::Uxs_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (hs>epsilon ? Uys*Uxs/hs : 0.); }

double
TG2_scheme::Uys_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ return (hs>epsilon ? Uys*Uxs/hs : 0.); }

double
TG2_scheme::Uys_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto vel_y = hs>epsilon ? Uys/hs : 0.;
  return (Uys*vel_y + grav*hs*hs/2. + grav*(1-r_coeff)*hs*hw/2.); 
}


// source terms
double
TG2_scheme::hw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto h = hw+hs;
  const auto Ux = Uxw + Uxs;
  const auto Uy = Uyw + Uys;
  const auto nw = h>epsilon ? hw/h : 0.;
  return (nw*erosion_coefficient*std::sqrt(Ux*Ux+Uy*Uy)); 
}

double
TG2_scheme::hs_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{ 
  const auto h = hw+hs;
  const auto Ux = Uxw + Uxs; 
  const auto Uy = Uyw + Uys;
  const auto ns = h>epsilon ? hs/h : 0.;
  return (ns*erosion_coefficient*std::sqrt(Ux*Ux+Uy*Uy)); 
}

double
TG2_scheme::Uxw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;

  const double vel_w_x = hw>epsilon ? Uxw/hw : 0.;
  const double vel_s_x = hs>epsilon ? Uxs/hs : 0.;

  const double C_d = hw>epsilon && h>epsilon ? n*ns/std::pow(n, m_coeff)/terminal_velocity*(density_s/density_w-1.)*grav : 0.;
  const double R_x = C_d*(vel_w_x-vel_s_x);


  return ( -h*R_x );
}

double
TG2_scheme::Uyw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;

  const double vel_w_y = hw>epsilon ? Uyw/hw : 0.;
  const double vel_s_y = hs>epsilon ? Uys/hs : 0.;

  const double C_d = (hw>epsilon && h>epsilon) ? n*ns/std::pow(n, m_coeff)/terminal_velocity*(density_s/density_w-1.)*grav : 0.;
  const double R_y = C_d*(vel_w_y-vel_s_y);

  return ( -h*R_y );
}

double
TG2_scheme::Uxs_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;
  const double bed_pressure = grav*h*ns*(1.-density_w/density_s); 
  const auto Ux = Uxw + Uxs; 
  const auto Uy = Uyw + Uys;
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::abs( vel_x );


  const double vel_w_x = hw>epsilon ? Uxw/hw : 0.;
  const double vel_s_x = hs>epsilon ? Uxs/hs : 0.;

  const double density = ns + n*density_w/density_s;

  const double C_d = (hw>epsilon && h>epsilon) ? n*ns/std::pow(n, m_coeff)/terminal_velocity*(1.-density_w/density_s)*grav : 0.;
  const double R_x = C_d*(vel_w_x-vel_s_x);

  //std::cout << h << " " << dhdx << " " << dZdx << " " << grav*h*(dZdx+dhdx) << std::endl;

  const double vel_x_sign = std::abs(Ux)>tolerance_sign ? Ux/std::abs(Ux) : Ux/tolerance_sign;
  //const double vel_x_sign = abs_vel>tolerance_sign ? vel_x/abs_vel : 0.;
  //const double vel_x_sign = abs_vel>tolerance_sign ? vel_x/abs_vel : vel_x/tolerance_sign;

  //const double bed_fric_contr = is_bed_friction ? vel_x_sign*(grav*abs_vel*abs_vel/turbulence_coeff + bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;

  const double bed_fric_contr_one = h*h>epsilon && is_bed_friction ? density*Ux*grav*std::abs(Ux)/turbulence_coeff/h/h : 0.; //is_bed_friction ? vel_x*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = is_bed_friction ? vel_x_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;

  return ( - bed_fric_contr_one - bed_fric_contr_two + h*R_x );
}

double
TG2_scheme::Uys_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys)
{
  const auto h = hw+hs;
  const auto n  = h>epsilon ? hw/h : 0.;
  const auto ns = h>epsilon ? hs/h : 0.;
  const double bed_pressure = grav*h*ns*(1.-density_w/density_s);
  const auto Ux = Uxw + Uxs;
  const auto Uy = Uyw + Uys;
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::abs( vel_y );



  const double vel_w_y = hw>epsilon ? Uyw/hw : 0.;
  const double vel_s_y = hs>epsilon ? Uys/hs : 0.;

  const double density = ns + n*density_w/density_s;

  const double C_d = (hw>epsilon && h>epsilon) ? n*ns/std::pow(n, m_coeff)/terminal_velocity*(1.-density_w/density_s)*grav : 0.;

  const double R_y = C_d*(vel_w_y-vel_s_y);

  const double vel_y_sign = std::abs(Uy)>tolerance_sign ? Uy/std::abs(Uy) : Uy/tolerance_sign;
  //const double vel_y_sign = abs_vel>tolerance_sign ? vel_y/abs_vel : vel_y/tolerance_sign;
  //const double vel_y_sign = (vel_y > ) ? 1.0 : (vel_y < 0) ? -1.0 : 0.0;

  //const double bed_fric_contr = is_bed_friction ? vel_y_sign*(grav*abs_vel*abs_vel/turbulence_coeff + bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;

  const double bed_fric_contr_one = ((h*h)>epsilon && is_bed_friction) ? density*Uy*grav*std::abs(Uy)/turbulence_coeff/h/h : 0.; //is_bed_friction ? vel_y*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = is_bed_friction ? vel_y_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;


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

