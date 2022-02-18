#include "Taylor_Galerkin_IMEX-RKC.h"
#include <algorithm>
#include <cassert>

TG2_scheme::TG2_scheme(Q1& sol,
                       Q1& sold,
                       Q1& soldd,
                       Q1& sold_rkc,
                       Q1& soldd_rkc,
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
                       const ordering& oUx,
                       const ordering& oUy,
                       const Q1& Z,
                       Q1& slope_x_node,
                       Q1& slope_y_node,
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
: sol(sol), sold(sold), soldd(soldd), sold_rkc(sold_rkc), soldd_rkc(soldd_rkc), incr(incr), incr_initial_source(incr_initial_source), incr_source(incr_source), incr_anti_diff(incr_anti_diff), stress_initial_step(stress_initial_step), stress_step(stress_step), P_plus(P_plus), P_minus(P_minus), spec_radius_nodal(spec_radius_nodal), sol_onehalf(sol_onehalf), mass(mass),
  ordh(oh), ordUx(oUx), ordUy(oUy), Z(Z), slope_x_node(slope_x_node), slope_y_node(slope_y_node), slope_x(slope_x), slope_y(slope_y), DELTAT(DELTAT), epsilon(h_min), is_non_reflBC(is_non_reflBC), is_bed_friction(is_bed_friction), is_stress_tensor(is_stress_tensor), grav(grav),
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
  h_cell_average    /= 4.;
  Ux_cell_average   /= 4.;
  Uy_cell_average   /= 4.;
  
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
  
  sol_onehalf[ordh    (index_quadrant)] = h_cell_average  - dt/2. * div_Fh_cell /area;
  sol_onehalf[ordUx   (index_quadrant)] = Ux_cell_average - dt/2. * div_FUx_cell/area + dt/2.*source_Ux_cell_average;
  sol_onehalf[ordUy   (index_quadrant)] = Uy_cell_average - dt/2. * div_FUy_cell/area + dt/2.*source_Uy_cell_average;
  


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

  double vel_rusanov_cell_x = 0., vel_rusanov_cell_y = 0.;
  for (int ii = 0; ii < 4; ++ii){
    const auto& hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    vel_rusanov_cell_x += hpoint>epsilon ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
    vel_rusanov_cell_y += hpoint>epsilon ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.; 

  grad_cell_h    = {.5 * ( (hdof [3] - hdof [2]) + (hdof [1] - hdof [0]) ), .5 * ( (hdof [2] - hdof [0]) + (hdof [3] - hdof [1]) )};
  grad_cell_Ux   = {.5 * ( (Uxdof[3] - Uxdof[2]) + (Uxdof[1] - Uxdof[0]) ), .5 * ( (Uxdof[2] - Uxdof[0]) + (Uxdof[3] - Uxdof[1]) )};
  grad_cell_Uy   = {.5 * ( (Uydof[3] - Uydof[2]) + (Uydof[1] - Uydof[0]) ), .5 * ( (Uydof[2] - Uydof[0]) + (Uydof[3] - Uydof[1]) )};

  

  const double & h_cell    = sol_onehalf[ordh    (index_quadrant)];
  const double & Ux_cell   = sol_onehalf[ordUx   (index_quadrant)];
  const double & Uy_cell   = sol_onehalf[ordUy   (index_quadrant)];


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
          break; // this just goes outside the jEdge cycle 
        }
      }

    }

    if (is_boundary_edge) // set boundary conditions
    { 
            
      auto h_cell_nei  = h_cell;
      auto Ux_cell_nei = Ux_cell;
      auto Uy_cell_nei = Uy_cell;

      Ux_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell + outward_normal_edge[1]*Uy_cell)*outward_normal_edge[0];
      Uy_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell + outward_normal_edge[1]*Uy_cell)*outward_normal_edge[1];

      const auto speed     = h_cell    >epsilon ? std::abs((Ux_cell    /h_cell    )*outward_normal_edge[0] + (Uy_cell    /h_cell    )*outward_normal_edge[1]) + std::sqrt(grav*h_cell    ) : 0.;
      const auto speed_nei = h_cell_nei>epsilon ? std::abs((Ux_cell_nei/h_cell_nei)*outward_normal_edge[0] + (Uy_cell_nei/h_cell_nei)*outward_normal_edge[1]) + std::sqrt(grav*h_cell_nei) : 0.;

      const auto smax = std::max(speed, speed_nei); 

      const auto flux_int_h  = .5*((h_flux_formula_x (h_cell, Ux_cell, Uy_cell)+h_flux_formula_x (h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (h_flux_formula_y (h_cell, Ux_cell, Uy_cell)+h_flux_formula_y (h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(h_cell_nei -h_cell );
      const auto flux_int_Ux = .5*((Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(Ux_cell_nei-Ux_cell);
      const auto flux_int_Uy = .5*((Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uy_cell_nei-Uy_cell);


      // occhio al segno qui!!!!!!!!!!!!!!!!, .5 is the base function evaluated in the middle, mid-point intergration
      if (! quadrant->is_hanging (i_1))
      {
        incr[ordh   (quadrant->gt (i_1))] += -edge_length*flux_int_h *.5;
        incr[ordUx  (quadrant->gt (i_1))] += -edge_length*flux_int_Ux*.5;
        incr[ordUy  (quadrant->gt (i_1))] += -edge_length*flux_int_Uy*.5;
      }
      else
      {
        incr [ordh  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_h *.5*.5;
        incr [ordh  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_h *.5*.5;
      
        incr [ordUx (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Ux*.5*.5;
        incr [ordUx (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Ux*.5*.5;
      
        incr [ordUy (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uy*.5*.5;
        incr [ordUy (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uy*.5*.5;
      }

      if (! quadrant->is_hanging (i_2))
      {
        incr[ordh   (quadrant->gt (i_2))] += -edge_length*flux_int_h *.5;
        incr[ordUx  (quadrant->gt (i_2))] += -edge_length*flux_int_Ux*.5;
        incr[ordUy  (quadrant->gt (i_2))] += -edge_length*flux_int_Uy*.5;
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
      }

    }

  }


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

    const auto h_    = der_coeffs_x[ii]*F_star_h_x +der_coeffs_y[ii]*F_star_h_y;
    const auto Ux_   = der_coeffs_x[ii]*F_star_Ux_x+der_coeffs_y[ii]*F_star_Ux_y;
    const auto Uy_   = der_coeffs_x[ii]*F_star_Uy_x+der_coeffs_y[ii]*F_star_Uy_y;

    const auto h_al  = der_coeffs_x[ii]*diff_term_h_x  + der_coeffs_y[ii]*diff_term_h_y; 
    const auto Ux_al = der_coeffs_x[ii]*diff_term_Ux_x + der_coeffs_y[ii]*diff_term_Ux_y;
    const auto Uy_al = der_coeffs_x[ii]*diff_term_Uy_x + der_coeffs_y[ii]*diff_term_Uy_y;

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
TG2_scheme::loop_step (const int& kk, const bool& isInitial) //
{

  const auto & h_c  = soldd_rkc.get_owned_data ()[kk  ];
  const auto & Ux_c = soldd_rkc.get_owned_data ()[kk+1];
  const auto & Uy_c = soldd_rkc.get_owned_data ()[kk+2]; 

    
  const auto h_s_  = 0.;
  const auto Ux_s_ = Ux_src_formula(h_c, Ux_c, Uy_c, slope_x_node.get_owned_data ()[int(kk/3)]); 
  const auto Uy_s_ = Uy_src_formula(h_c, Ux_c, Uy_c, slope_y_node.get_owned_data ()[int(kk/3)]); 

  if (isInitial)
  {
    incr_initial_source.get_owned_data ()[kk  ] = h_s_; 
    incr_initial_source.get_owned_data ()[kk+1] = Ux_s_;
    incr_initial_source.get_owned_data ()[kk+2] = Uy_s_;    
  }
  else
  {
    incr_source.get_owned_data ()[kk  ] = h_s_; 
    incr_source.get_owned_data ()[kk+1] = Ux_s_;
    incr_source.get_owned_data ()[kk+2] = Uy_s_;
  }
  


  /*
  // look at tmesh.h
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 
  
  Dx = quadrant->p(0,1) - quadrant->p(0,0);
  Dy = quadrant->p(1,2) - quadrant->p(1,0);
  area = Dx * Dy;

  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){
      hdof[ii]    = soldd_rkc [ordh    (quadrant->gt (ii))];
      Uxdof[ii]   = soldd_rkc [ordUx   (quadrant->gt (ii))];
      Uydof[ii]   = soldd_rkc [ordUy   (quadrant->gt (ii))];
      
      isdof_or_hanging[ii] = 1.;

    } else {
      hdof[ii]    = .5 * (soldd_rkc [ordh  (quadrant->gparent(0,ii))] +
                          soldd_rkc [ordh  (quadrant->gparent(1,ii))]);
      Uxdof[ii]   = .5 * (soldd_rkc [ordUx (quadrant->gparent(0,ii))] +
                          soldd_rkc [ordUx (quadrant->gparent(1,ii))]);
      Uydof[ii]   = .5 * (soldd_rkc [ordUy (quadrant->gparent(0,ii))] +
                          soldd_rkc [ordUy (quadrant->gparent(1,ii))]);
      
      isdof_or_hanging[ii] = .5;
    }
  }

  for (int ii = 0; ii < 4; ++ii){
    
    const auto h_s_  = 0.;
    const auto Ux_s_ = .25*area*isdof_or_hanging[ii]*Ux_src_formula(hdof[ii], Uxdof[ii], Uydof[ii], slope_x[index_quadrant]); 
    const auto Uy_s_ = .25*area*isdof_or_hanging[ii]*Uy_src_formula(hdof[ii], Uxdof[ii], Uydof[ii], slope_y[index_quadrant]); 

    if (!isInitial)
    {
      if (! quadrant->is_hanging (ii)){

        incr_source [ordh  (quadrant->gt (ii))] += h_s_;
        incr_source [ordUx (quadrant->gt (ii))] += Ux_s_;
        incr_source [ordUy (quadrant->gt (ii))] += Uy_s_;

      } else {

        incr_source [ordh  (quadrant->gparent(0,ii))] += h_s_;
        incr_source [ordh  (quadrant->gparent(1,ii))] += h_s_;

        incr_source [ordUx (quadrant->gparent(0,ii))] += Ux_s_;
        incr_source [ordUx (quadrant->gparent(1,ii))] += Ux_s_;

        incr_source [ordUy (quadrant->gparent(0,ii))] += Uy_s_;
        incr_source [ordUy (quadrant->gparent(1,ii))] += Uy_s_;

      }
    }
    else
    {
      if (! quadrant->is_hanging (ii)){

        incr_initial_source [ordh  (quadrant->gt (ii))] += h_s_;
        incr_initial_source [ordUx (quadrant->gt (ii))] += Ux_s_;
        incr_initial_source [ordUy (quadrant->gt (ii))] += Uy_s_;

      } else {

        incr_initial_source [ordh  (quadrant->gparent(0,ii))] += h_s_;
        incr_initial_source [ordh  (quadrant->gparent(1,ii))] += h_s_;

        incr_initial_source [ordUx (quadrant->gparent(0,ii))] += Ux_s_;
        incr_initial_source [ordUx (quadrant->gparent(1,ii))] += Ux_s_;

        incr_initial_source [ordUy (quadrant->gparent(0,ii))] += Uy_s_;
        incr_initial_source [ordUy (quadrant->gparent(1,ii))] += Uy_s_;

      }
    }
  }*/


}


void
TG2_scheme::compute_stress_slope (tmesh::quadrant_iterator quadrant, const bool& isInitial)
{
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 
  const std::array<double,4> zeros_array = {0, 0, 0, 0};
  std::array<std::array<double,4>, 4 > delta_incr_mat = {zeros_array, zeros_array, zeros_array, zeros_array};

  std::array<int,4> bimpp_to_rev_ord = {0, 1, 3, 2};

  Dx = quadrant->p(0,1) - quadrant->p(0,0);
  Dy = quadrant->p(1,2) - quadrant->p(1,0);
  area = Dx*Dy;

  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){
      hdof[ii]    = sold_rkc [ordh    (quadrant->gt (ii))];
      Uxdof[ii]   = sold_rkc [ordUx   (quadrant->gt (ii))];
      Uydof[ii]   = sold_rkc [ordUy   (quadrant->gt (ii))];
      
      isdof_or_hanging[ii] = 1.;

      delta_incr_mat[ii][ii] = tol_incr;
    } else {
      hdof[ii]    = .5 * (sold_rkc [ordh  (quadrant->gparent(0,ii))] +
                          sold_rkc [ordh  (quadrant->gparent(1,ii))]);
      Uxdof[ii]   = .5 * (sold_rkc [ordUx (quadrant->gparent(0,ii))] +
                          sold_rkc [ordUx (quadrant->gparent(1,ii))]);
      Uydof[ii]   = .5 * (sold_rkc [ordUy (quadrant->gparent(0,ii))] +
                          sold_rkc [ordUy (quadrant->gparent(1,ii))]);
      
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

  //std::cout << sigma_stress[0] << " " << sigma_stress[1] << " " << sigma_stress[2] << " " << sigma_stress.size() << std::endl;

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


  /*
  // boundary conditions, here are for the diffusion term only, if commented means null diffusive fluxes
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
          break; // this just goes outside the jEdge cycle 
        }
      }

    }

    
    if (is_boundary_edge) // set boundary conditions
    { 
      const double middle_D_U = (D_U[i_1]+D_U[i_2])*.5; 

      const double Integral_b = 1./6*(D_U[i_1]+4.*(middle_D_U*.5)+D_U[i_2]);

      const auto flux_int_h  = 0.;
      const auto flux_int_Ux = (sigma_stress[0]*outward_normal_edge[0] + sigma_stress[2]*outward_normal_edge[1])*Integral_b;
      const auto flux_int_Uy = (sigma_stress[2]*outward_normal_edge[0] + sigma_stress[1]*outward_normal_edge[1])*Integral_b;

      // imponiamo il flusso degli sforzi nullo al bordo
      // .5 is the base function evaluated in the middle, mid-point intergration
      if (isInitial)
      {
        if (! quadrant->is_hanging (i_1))
        {
          stress_initial_step[ordh   (quadrant->gt (i_1))] += -edge_length*flux_int_h ;
          stress_initial_step[ordUx  (quadrant->gt (i_1))] += -edge_length*flux_int_Ux;
          stress_initial_step[ordUy  (quadrant->gt (i_1))] += -edge_length*flux_int_Uy;
        }
        else
        {
          stress_initial_step [ordh  (quadrant->gparent(0,i_1))] += -edge_length*flux_int_h *.5;
          stress_initial_step [ordh  (quadrant->gparent(1,i_1))] += -edge_length*flux_int_h *.5;

          stress_initial_step [ordUx (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Ux*.5;
          stress_initial_step [ordUx (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Ux*.5;

          stress_initial_step [ordUy (quadrant->gparent(0,i_1))] += -edge_length*flux_int_Uy*.5;
          stress_initial_step [ordUy (quadrant->gparent(1,i_1))] += -edge_length*flux_int_Uy*.5;
        }
      }
      else
      {
        if (! quadrant->is_hanging (i_2))
        {
          stress_step[ordh   (quadrant->gt (i_2))] += -edge_length*flux_int_h ;
          stress_step[ordUx  (quadrant->gt (i_2))] += -edge_length*flux_int_Ux;
          stress_step[ordUy  (quadrant->gt (i_2))] += -edge_length*flux_int_Uy;
        }
        else
        {
          // il secondo .5 è per hanging nodes
          stress_step [ordh  (quadrant->gparent(0,i_2))] += -edge_length*flux_int_h *.5;
          stress_step [ordh  (quadrant->gparent(1,i_2))] += -edge_length*flux_int_h *.5;

          stress_step [ordUx (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Ux*.5;
          stress_step [ordUx (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Ux*.5;

          stress_step [ordUy (quadrant->gparent(0,i_2))] += -edge_length*flux_int_Uy*.5;
          stress_step [ordUy (quadrant->gparent(1,i_2))] += -edge_length*flux_int_Uy*.5;
        }
      }

    }

  }
  */


  for (int ii = 0; ii < 4; ++ii){

    const double den1 = ii<2    ? 2. : 1.;
    const double den2 = ii<2    ? 1. : 2.;
    const double den3 = ii%2==1 ? 2. : 1.;
    const double den4 = ii%2==1 ? 1. : 2.;


    const double contribution_exact = (D_U[0]/den2+D_U[1]/den2+D_U[2]/den1+D_U[3]/den1);    

    const double h_  = 0.;
    const double Ux_ = der_coeffs_x[ii]*(1./3.)*sigma_stress[0]*contribution_exact + der_coeffs_y[ii]*(1./3.)*sigma_stress[2]*contribution_exact;
    const double Uy_ = der_coeffs_x[ii]*(1./3.)*sigma_stress[2]*contribution_exact + der_coeffs_y[ii]*(1./3.)*sigma_stress[1]*contribution_exact;

    double h_s = 0., Ux_s = 0., Uy_s = 0.;
    for (int jj = 0; jj < 4; ++jj)
    {
      //std::cout << sigma_stress_incr_Ux[jj][0] << " " << sigma_stress[0] << std::endl;
      //std::cout << der_coeffs_x[ii]*(1./3.)*(sigma_stress_incr_Ux[jj][0]-sigma_stress[0])/tol_incr*contribution_exact << std::endl;

      Ux_s += std::abs( der_coeffs_x[ii]*(1./3.)*(sigma_stress_incr_Ux[jj][0]-sigma_stress[0])/tol_incr*contribution_exact + der_coeffs_y[ii]*(1./3.)*(sigma_stress_incr_Ux[jj][2]-sigma_stress[2])/tol_incr*contribution_exact );
      Ux_s += std::abs( der_coeffs_x[ii]*(1./3.)*(sigma_stress_incr_Uy[jj][0]-sigma_stress[0])/tol_incr*contribution_exact + der_coeffs_y[ii]*(1./3.)*(sigma_stress_incr_Uy[jj][2]-sigma_stress[2])/tol_incr*contribution_exact );

      Uy_s += std::abs( der_coeffs_x[ii]*(1./3.)*(sigma_stress_incr_Ux[jj][2]-sigma_stress[2])/tol_incr*contribution_exact + der_coeffs_y[ii]*(1./3.)*(sigma_stress_incr_Ux[jj][1]-sigma_stress[1])/tol_incr*contribution_exact );
      Uy_s += std::abs( der_coeffs_x[ii]*(1./3.)*(sigma_stress_incr_Uy[jj][2]-sigma_stress[2])/tol_incr*contribution_exact + der_coeffs_y[ii]*(1./3.)*(sigma_stress_incr_Uy[jj][1]-sigma_stress[1])/tol_incr*contribution_exact );
    }
    

    //std::cout << sigma_stress[0] << " " << sigma_stress[1] << " " << sigma_stress[2] << " " << Ux_ << " " << Uy_ << std::endl;

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

        /*
        spec_radius_nodal [ordh   (quadrant->gt (ii))] += h_s;
        spec_radius_nodal [ordUx  (quadrant->gt (ii))] += Ux_s;
        spec_radius_nodal [ordUy  (quadrant->gt (ii))] += Uy_s;*/

      } else {

        stress_step [ordh  (quadrant->gparent(0,ii))] += h_;
        stress_step [ordh  (quadrant->gparent(1,ii))] += h_;

        stress_step [ordUx (quadrant->gparent(0,ii))] += Ux_;
        stress_step [ordUx (quadrant->gparent(1,ii))] += Ux_;

        stress_step [ordUy (quadrant->gparent(0,ii))] += Uy_;
        stress_step [ordUy (quadrant->gparent(1,ii))] += Uy_;

        /*
        spec_radius_nodal [ordh  (quadrant->gparent(0,ii))] += h_s;
        spec_radius_nodal [ordh  (quadrant->gparent(1,ii))] += h_s;

        spec_radius_nodal [ordUx (quadrant->gparent(0,ii))] += Ux_s;
        spec_radius_nodal [ordUx (quadrant->gparent(1,ii))] += Ux_s;

        spec_radius_nodal [ordUy (quadrant->gparent(0,ii))] += Uy_s;
        spec_radius_nodal [ordUy (quadrant->gparent(1,ii))] += Uy_s;*/

      }
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
TG2_scheme::rkc(const int& j, const int& s, const int& kk)
{
  // kk is the current owned node

  double v_x, v_y;

  const auto & S_x  = slope_x_node.get_owned_data ()[int(kk/3)];
  const auto & S_y  = slope_y_node.get_owned_data ()[int(kk/3)]; 

  int count; 

  const double tolerance = 1.e-4;
  const int Nmax = 1e3;

  double error; 

  if (j == 1)
  {
    //std::cout << stress_initial_step.get_owned_data ()[kk+1] << std::endl;

    v_x = sol.get_owned_data ()[kk+1] + mu_fun_tilde(1, s)*dt*(incr.get_owned_data ()[kk+1] + stress_initial_step.get_owned_data ()[kk+1])/mass.get_owned_data ()[kk+1];
    v_y = sol.get_owned_data ()[kk+2] + mu_fun_tilde(1, s)*dt*(incr.get_owned_data ()[kk+2] + stress_initial_step.get_owned_data ()[kk+2])/mass.get_owned_data ()[kk+2]; 
  }
  else
  {
    //std::cout << "aa" << std::endl;

    v_x = (1. - mu_fun(j, s) - v_fun(j, s))*sold.get_owned_data ()[kk+1] + mu_fun(j, s)*sold_rkc.get_owned_data ()[kk+1] + 
    v_fun(j, s)*soldd_rkc.get_owned_data ()[kk+1] + mu_fun_tilde(j, s)*dt*(incr.get_owned_data ()[kk+1] + stress_step.get_owned_data ()[kk+1])/mass.get_owned_data ()[kk+1] + 
    gamma_tilde_fun(j, s)*dt*(incr.get_owned_data ()[kk+1] + stress_initial_step.get_owned_data ()[kk+1])/mass.get_owned_data ()[kk+1] + 
    (gamma_tilde_fun(j, s) - (1. - mu_fun(j, s) - v_fun(j, s))*mu_fun_tilde(1, s))*dt*incr_initial_source.get_owned_data ()[kk+1] - v_fun(j, s)*mu_fun_tilde(1, s)*dt*incr_source.get_owned_data ()[kk+1];

    v_y = (1. - mu_fun(j, s) - v_fun(j, s))*sold.get_owned_data ()[kk+2] + mu_fun(j, s)*sold_rkc.get_owned_data ()[kk+2] + 
    v_fun(j, s)*soldd_rkc.get_owned_data ()[kk+2] + mu_fun_tilde(j, s)*dt*(incr.get_owned_data ()[kk+2] + stress_step.get_owned_data ()[kk+2])/mass.get_owned_data ()[kk+2] + 
    gamma_tilde_fun(j, s)*dt*(incr.get_owned_data ()[kk+2] + stress_initial_step.get_owned_data ()[kk+2])/mass.get_owned_data ()[kk+2] +
    (gamma_tilde_fun(j, s) - (1. - mu_fun(j, s) - v_fun(j, s))*mu_fun_tilde(1, s))*dt*incr_initial_source.get_owned_data ()[kk+2] - v_fun(j, s)*mu_fun_tilde(1, s)*dt*incr_source.get_owned_data ()[kk+2];
  }


  // Ux
  count = -1;
  error = tolerance + 1;
  while (count++<Nmax && error>tolerance)
  {
    const auto & h_c  = sol.get_owned_data ()[kk  ];
    const auto & Ux_c = sol.get_owned_data ()[kk+1];
    const auto & Uy_c = sol.get_owned_data ()[kk+2];

    const auto delta_Ux = (- Ux_c + v_x + mu_fun_tilde(1, s)*dt*Ux_src_formula(h_c, Ux_c, Uy_c, S_x))/Ux_jac_source(h_c, Ux_c, Uy_c, s);

    error = std::abs(delta_Ux);

    //std::cout << v_x << " " << delta_Ux << " " << sol.get_owned_data ()[kk+1] << " " << Ux_src_formula(h_c, Ux_c, Uy_c, S_x) << " " << S_x << " " << is_bed_friction << std::endl;

    sol.get_owned_data ()[kk+1] += delta_Ux;

    //if (count==1) std::cout << count << " " << delta_Ux << std::endl; 

  }


  // Uy
  count = -1;
  error = tolerance + 1;
  while (count++<Nmax && error>tolerance)
  {
    const auto & h_c  = sol.get_owned_data ()[kk  ];
    const auto & Ux_c = sol.get_owned_data ()[kk+1];
    const auto & Uy_c = sol.get_owned_data ()[kk+2];

    const auto delta_Uy = (- Uy_c + v_y + mu_fun_tilde(1, s)*dt*Uy_src_formula(h_c, Ux_c, Uy_c, S_y))/Uy_jac_source(h_c, Ux_c, Uy_c, s);

    error = std::abs(delta_Uy);

    sol.get_owned_data ()[kk+2] += delta_Uy;

  }

}


double
TG2_scheme::Ux_jac_source(const double& h, const double& Ux, const double& Uy, const int& s)
{
  return((h>epsilon && is_bed_friction) ? 1.+mu_fun_tilde(1, s)*dt*grav/turbulence_coeff/h/h*2.*std::abs(Ux) : 1.);
}

double
TG2_scheme::Uy_jac_source(const double& h, const double& Ux, const double& Uy, const int& s)
{
  return((h>epsilon && is_bed_friction) ? 1.+mu_fun_tilde(1, s)*dt*grav/turbulence_coeff/h/h*2.*std::abs(Uy) : 1.);
}

/*
double
TG2_scheme::c_fun(const int& j, const int& s)
{
  double c;

  if (j == 0)
  {
    c = 0;
  }
  else if (j == 1)
  {
    c = (j*j + 2*j)/(s*s - 1)/w_fun_0(s);
  }
  else if (j == s)
  {
    c = 1.;
  }
  else
  {
    c = (j*j - 1)/(s*s - 1);
  }

  return(c);
}
*/

double
TG2_scheme::mu_fun(const int& j, const int& s)
{
  const auto mu = 2.*b_fun(j, s)*w_fun_0(s)/b_fun(j-1, s);
  return(mu);
}

double
TG2_scheme::v_fun(const int& j, const int& s)
{
  const auto v = -b_fun(j, s)/b_fun(j-2, s);
  return(v);
}

double
TG2_scheme::gamma_tilde_fun(const int& j, const int& s)
{
  const auto w0 = w_fun_0(s);
  const auto gamma_tilde = -(1-b_fun(j-1, s)*T_fun(j-1, w0))*mu_fun_tilde(j,s);

  return(gamma_tilde);
}

double
TG2_scheme::T_fun_second(const int& s, const double& x)
{
  double T;

  if (s == 0)
  {
    T = 0;
  }
  else if (s == 1)
  {
    T = 0;
  }
  else
  {
    T = 2.*T_fun_prime(s-1,x) + 2.*T_fun_prime(s-1,x) + 2.*x*T_fun_second(s-1,x) - T_fun_second(s-2,x);
  }

  return(T);
}

double
TG2_scheme::T_fun_prime(const int& s, const double& x)
{
  double T;

  if (s == 0)
  {
    T = 0;
  }
  else if (s == 1)
  {
    T = 1;
  }
  else
  {
    T = 2.*T_fun(s-1,x) + 2.*x*T_fun_prime(s-1,x) - T_fun_prime(s-2,x);
  }

  return(T);
}

double
TG2_scheme::T_fun(const int& s, const double& x)
{
  double T;

  if (s == 0)
  {
    T = 1;
  }
  else if (s == 1)
  {
    T = x;
  }
  else
  {
    T = 2.*x*T_fun(s-1,x) - T_fun(s-2,x);
  }

  return(T);
}

double
TG2_scheme::w_fun_0(const int& s)
{
  const double epsilon_w = 2./13.;
  const double w = 1. + epsilon_w/s/s;
  return(w);
}

double
TG2_scheme::w_fun_1(const int& s)
{
  const auto w0 = w_fun_0(s);
  const auto w = T_fun_prime(s,w0)/T_fun_second(s,w0);

  return(w);
}

double
TG2_scheme::b_fun(const int& j, const int& s)
{
  const auto w0 = w_fun_0(s);

  double b;

  if (j <= 1)
  {
    b = T_fun_second(2,w0)/T_fun_prime(2,w0)/T_fun_prime(2,w0);
  }
  else
  {
    b = T_fun_second(j,w0)/T_fun_prime(j,w0)/T_fun_prime(j,w0);
  }

  return(b);


}



double
TG2_scheme::mu_fun_tilde (const int& j, const int& s)
{
  double mu_tilde;

  if (j == 1)
  {
    mu_tilde = b_fun(1, s)*w_fun_1(s);
  }
  else
  {
    mu_tilde = 2.*b_fun(j, s)*w_fun_1(s)/b_fun(j-1, s);
  }

  return(mu_tilde);
    
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

  double viscos = second_invariant!=0 ? 2.*fluid_viscosity + yield_shear_stress/kinetic_energergy_associated*(1. - std::exp(-regularization_parameter*kinetic_energergy_associated)) : 0.;
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
    std::cout << "Two valid roots, look at compute_nodal_def_grad function, " << zeta_1 << " " << zeta_2 << ", STOP!" << std::endl;
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
TG2_scheme::h_src_formula (const double& h, const double& Ux, const double& Uy)
{ return (0.); }

double
TG2_scheme::Ux_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdx)
{
  const double bed_pressure = grav*h + surface_pressure/density; 
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::sqrt( vel_x*vel_x + vel_y*vel_y );

  const double vel_x_sign = abs_vel!=0 ? vel_x/abs_vel : 0.;

  //const double bed_fric_contr = is_bed_friction ? vel_x_sign*(grav*abs_vel*abs_vel/turbulence_coeff + bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;

  const double bed_fric_contr_one = is_bed_friction ? vel_x*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = is_bed_friction ? vel_x_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;


  return (- grav*h*dZdx - bed_fric_contr_one - bed_fric_contr_two);
}

double
TG2_scheme::Uy_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdy)
{
  const double bed_pressure = grav*h + surface_pressure/density;
  const double vel_x = h>epsilon ? Ux/h : 0.;
  const double vel_y = h>epsilon ? Uy/h : 0.;
  const double abs_vel = std::sqrt( vel_x*vel_x + vel_y*vel_y );

  const double vel_y_sign = abs_vel!=0 ? vel_y/abs_vel : 0.;

  //const double bed_fric_contr = is_bed_friction ? vel_y_sign*(grav*abs_vel*abs_vel/turbulence_coeff + bed_pressure*std::tan(bed_friction_angle_rad)) : 0.;

  const double bed_fric_contr_one = is_bed_friction ? vel_y*grav*abs_vel/turbulence_coeff : 0.;
  const double bed_fric_contr_two = is_bed_friction ? vel_y_sign*bed_pressure*std::tan(bed_friction_angle_rad) : 0.;

  //std::cout << bed_pressure << std::endl;

  return (- grav*h*dZdy - bed_fric_contr_one - bed_fric_contr_two);
}



