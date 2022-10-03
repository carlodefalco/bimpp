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
		       const bool& is_max_time_step_from_CFL_with_diffusion,
                       const double& grav,
                       const double& density,
                       const double& turbulence_coeff,
                       const double& surface_pressure, 
                       const double& bed_friction_angle_rad,
                       const double& fluid_viscosity,
                       const double& yield_shear_stress)
: sol(sol), sold(sold), soldd(soldd), incr(incr), incr_anti_diff(incr_anti_diff), P_plus(P_plus), P_minus(P_minus), sol_onehalf(sol_onehalf), mass(mass),
  ordh(oh), ordUx(oUx), ordUy(oUy), Z(Z), slope_x(slope_x), slope_y(slope_y), DELTAT(DELTAT), epsilon(h_min), is_non_reflBC(is_non_reflBC), is_bed_friction(is_bed_friction), is_stress_tensor(is_stress_tensor), is_max_time_step_from_CFL_with_diffusion(is_max_time_step_from_CFL_with_diffusion), grav(grav),
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
    
    const auto vel_rusanov_cell_x = hpoint>epsilon ? std::max(std::abs(Uxdof[ii]/hpoint)+celerity, is_stress_tensor*is_max_time_step_from_CFL_with_diffusion ? 2*fluid_viscosity/Dx : 0.) : 0.;
    const auto vel_rusanov_cell_y = hpoint>epsilon ? std::max(std::abs(Uydof[ii]/hpoint)+celerity, is_stress_tensor*is_max_time_step_from_CFL_with_diffusion ? 2*fluid_viscosity/Dy : 0.) : 0.;
    
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
  nu_htot += Nu_hmean_cell*(time-timed)*(time-timed); // the dimension is L^2, this is eta^2
  
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
  
  sol_onehalf[ordh    (index_quadrant)] = h_cell_average  - dt/2.*div_Fh_cell /area;
  sol_onehalf[ordUx   (index_quadrant)] = Ux_cell_average - dt/2.*div_FUx_cell/area + dt/2.*source_Ux_cell_average;
  sol_onehalf[ordUy   (index_quadrant)] = Uy_cell_average - dt/2.*div_FUy_cell/area + dt/2.*source_Uy_cell_average;
  


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

  grad_cell_ux   = {.5 * ( (hdof[3]>epsilon ? Uxdof[3]/hdof[3] : 0. - hdof[2]>epsilon ? Uxdof[2]/hdof[2] : 0.) + (hdof[1]>epsilon ? Uxdof[1]/hdof[1] : 0. - hdof[0]>epsilon ? Uxdof[0]/hdof[0] : 0.) )/Dx, .5 * ( (hdof[2]>epsilon ? Uxdof[2]/hdof[2] : 0. - hdof[0]>epsilon ? Uxdof[0]/hdof[0] : 0.) + (hdof[3]>epsilon ? Uxdof[3]/hdof[3] : 0. - hdof[1]>epsilon ? Uxdof[1]/hdof[1] : 0.) )/Dy};
  grad_cell_uy   = {.5 * ( (hdof[3]>epsilon ? Uydof[3]/hdof[3] : 0. - hdof[2]>epsilon ? Uydof[2]/hdof[2] : 0.) + (hdof[1]>epsilon ? Uydof[1]/hdof[1] : 0. - hdof[0]>epsilon ? Uydof[0]/hdof[0] : 0.) )/Dx, .5 * ( (hdof[2]>epsilon ? Uydof[2]/hdof[2] : 0. - hdof[0]>epsilon ? Uydof[0]/hdof[0] : 0.) + (hdof[3]>epsilon ? Uydof[3]/hdof[3] : 0. - hdof[1]>epsilon ? Uydof[1]/hdof[1] : 0.) )/Dy};


  

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





  const auto diff_term_h_x  = grad_cell_h [0]*vel_rusanov_cell_y;
  const auto diff_term_h_y  = grad_cell_h [1]*vel_rusanov_cell_x;

  const auto diff_term_Ux_x = grad_cell_Ux[0]*vel_rusanov_cell_y;
  const auto diff_term_Ux_y = grad_cell_Ux[1]*vel_rusanov_cell_x;

  const auto diff_term_Uy_x = grad_cell_Uy[0]*vel_rusanov_cell_y;
  const auto diff_term_Uy_y = grad_cell_Uy[1]*vel_rusanov_cell_x;
  


  const auto F_star_h_x  = h_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_h_x;
  const auto F_star_h_y  = h_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_h_y;

  const auto F_star_Ux_x = Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_x;
  const auto F_star_Ux_y = Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_y;

  const auto F_star_Uy_x = Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_x;
  const auto F_star_Uy_y = Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_y;





  // compute the cell sigma_stress
  sigma_stress = is_stress_tensor ? compute_cell_stress (Uxdof[0], Uxdof[1], Uxdof[2], Uxdof[3], Uydof[0], Uydof[1], Uydof[2], Uydof[3]) : sigma_stress;

  for (int ii = 0; ii < 4; ++ii){
    D_U[ii] = U_stress_formula(hdof[ii], Uxdof[ii], Uydof[ii]);
  }

  for (int ii = 0; ii < 4; ++ii){ 

    const double den1 = ii<2    ? 2. : 1.;
    const double den2 = ii<2    ? 1. : 2.;
    const double den3 = ii%2==1 ? 2. : 1.;
    const double den4 = ii%2==1 ? 1. : 2.;

    const double contribution_exact = (D_U[0]/den2+D_U[1]/den2+D_U[2]/den1+D_U[3]/den1);    

    const auto h_  = der_coeffs_x[ii]*F_star_h_x +der_coeffs_y[ii]*F_star_h_y;

    const auto Ux_ = der_coeffs_x[ii]*F_star_Ux_x+der_coeffs_y[ii]*F_star_Ux_y + 
                     der_coeffs_x[ii]*(1./3.)*sigma_stress[0]*contribution_exact + der_coeffs_y[ii]*(1./3.)*sigma_stress[2]*contribution_exact +
                     .25*area*isdof_or_hanging[ii]*Ux_src_formula(h_cell, Ux_cell, Uy_cell, slope_x[index_quadrant]);

    const auto Uy_ = der_coeffs_x[ii]*F_star_Uy_x+der_coeffs_y[ii]*F_star_Uy_y +  
                     der_coeffs_x[ii]*(1./3.)*sigma_stress[2]*contribution_exact + der_coeffs_y[ii]*(1./3.)*sigma_stress[1]*contribution_exact + 
                     .25*area*isdof_or_hanging[ii]*Uy_src_formula(h_cell, Ux_cell, Uy_cell, slope_y[index_quadrant]);
    
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


// // stress functions
// double
// TG2_scheme::Ux_stress_formula_x (const double& h, const double& Ux, const double& Uy)
// { return (-sigma_stress[0]*h/density); }

// double
// TG2_scheme::Ux_stress_formula_y (const double& h, const double& Ux, const double& Uy)
// { return (-sigma_stress[2]*h/density); }

// double
// TG2_scheme::Uy_stress_formula_x (const double& h, const double& Ux, const double& Uy)
// { return (-sigma_stress[2]*h/density); }

// double
// TG2_scheme::Uy_stress_formula_y (const double& h, const double& Ux, const double& Uy)
// { return (-sigma_stress[1]*h/density); }




// std::array<double,3>
// TG2_scheme::compute_nodal_stress (const double& h, const double& Ux, const double& Uy, const std::array<double,2>& grad_cell_ux, const std::array<double,2>& grad_cell_uy)
// {
//   // compute \sigma_xx, ...

//   // def_grad = [D11, D22, D33, D12, D23, D31]
//   // sigma = [sigma_11, sigma_22, sigma_12]

//   std::array<double,6> def_grad = compute_nodal_def_grad (h, Ux, Uy, grad_cell_ux, grad_cell_uy);

//   double second_invariant = 0.;
//   for (int i_def = 0; i_def < 6; i_def++)
//   {
//     second_invariant += def_grad[i_def]*def_grad[i_def];
//   }
//   second_invariant *= .5;

//   //const double viscos = second_invariant!=0 ? yield_shear_stress/std::sqrt(second_invariant) + 2*fluid_viscosity : 0.;

//   const double kinetic_energergy_associated = std::sqrt(second_invariant);

//   double viscos = second_invariant!=0 ? 2.*fluid_viscosity + yield_shear_stress/kinetic_energergy_associated*(1. - std::exp(-regularization_parameter*kinetic_energergy_associated)) : 0.;
//   viscos = yield_shear_stress==0 ? 2.*fluid_viscosity : viscos;

// //std::cout << viscos << std::endl;//def_grad[0] << " " << def_grad[1] << " " << def_grad[2] << std::endl;
//   return(std::array<double,3>{{viscos*def_grad[0], viscos*def_grad[1], viscos*def_grad[3]}});
// }

// std::array<double,6>
// TG2_scheme::compute_nodal_def_grad (const double& h, const double& Ux, const double& Uy, const std::array<double,2>& grad_cell_ux, const std::array<double,2>& grad_cell_uy)
// {

//   // def_grad = [D11, D22, D33, D12, D23, D31]

//   // compute \zeta
//   const double vel_x = h>epsilon ? Ux/h : 0.;
//   const double vel_y = h>epsilon ? Uy/h : 0.;
//   const double abs_vel = std::sqrt( vel_x*vel_x + vel_y*vel_y );
//   const double aa = h>epsilon ? 6*fluid_viscosity*abs_vel/h/yield_shear_stress : 0.;

//   const double a = 3./2.;
//   const double c = 65./32.;
//   const double b = -(114./32.+aa);
//   const double Delta = b*b-4*a*c;

//   const double zeta_1 = (-b + std::sqrt(Delta))/2./a;
//   const double zeta_2 = (-b - std::sqrt(Delta))/2./a;

//   if ( std::abs(zeta_1 - .5)<=.5 && std::abs(zeta_2 - .5)<=.5)
//   {
//     std::cout << "Two valid roots, look at compute_nodal_def_grad, " << zeta_1 << " " << zeta_2 << ", STOP!" << std::endl;
//     exit(1.);
//   }

//   const double zeta = std::abs(zeta_1 - .5)<=.5 ? zeta_1 : zeta_2;

//   const auto & partial_x_ux = grad_cell_ux[0];
//   const auto & partial_y_ux = grad_cell_ux[1];
//   const auto   partial_z_ux = h>epsilon ? 3./(2.+zeta)*vel_x/h : 0.;

//   const auto & partial_x_uy = grad_cell_uy[0];
//   const auto & partial_y_uy = grad_cell_uy[1];
//   const auto   partial_z_uy = h>epsilon ? 3./(2.+zeta)*vel_y/h : 0.;

//   const auto partial_x_uz = 0.; // steady state simple shear flow
//   const auto partial_y_uz = 0.; // steady state simple shear flow
//   const auto partial_z_uz = -(grad_cell_ux[0]+grad_cell_uy[1]);
  
//   return(std::array<double,6>{{partial_x_ux, partial_y_uy, partial_z_uz, .5*(partial_x_uy+partial_y_ux), .5*(partial_z_uy+partial_y_uz), .5*(partial_z_ux+partial_x_uz)}});
// }


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

  //if (abs_vel>100)//(std::abs(bed_fric_contr_one)>std::abs(bed_fric_contr_two))
  //{
    //std::cout << abs_vel << " " << Ux << " " << bed_fric_contr_one << " " << vel_x*grav/turbulence_coeff << " " << bed_fric_contr_two << " " << bed_fric_contr_one+bed_fric_contr_two << std::endl;
    //exit(1);
  //}

  return (-grav*h*dZdx - bed_fric_contr_one - bed_fric_contr_two);
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

  return (-grav*h*dZdy - bed_fric_contr_one - bed_fric_contr_two);
}



