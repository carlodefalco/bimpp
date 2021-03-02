#include "Taylor_Galerkin.h"



TG2_scheme::TG2_scheme(const Q1& state,
                       Q0& state_onehalf,
                       const ordering& oh,
                       const ordering& oUx,
                       const ordering& oUy,
                       const Q1& elevation,
                       const double& DELTAT)
: state_vector(state), sol_onehalf(state_onehalf), ordh(oh), ordUx(oUx), ordUy(oUy), Z(elevation), DELTAT(DELTAT)
{ gn_elements = int(std::round(sol_onehalf.size()/3.)); local_extrema.resize(3*gn_elements);}



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
      hdof[ii]  = state_vector [ordh (quadrant->gt (ii) )];
      Uxdof[ii] = state_vector [ordUx (quadrant->gt (ii) )];
      Uydof[ii] = state_vector [ordUy (quadrant->gt (ii) )];
      
    }
    else
    {
      hdof[ii]  = .5 * (state_vector [ordh (quadrant->gparent (0, ii) )] +
                        state_vector [ordh (quadrant->gparent (1, ii) )]);
      Uxdof[ii] = .5 * (state_vector [ordUx (quadrant->gparent (0, ii) )] +
                        state_vector [ordUx (quadrant->gparent (1, ii) )]);
      Uydof[ii] = .5 * (state_vector [ordUy (quadrant->gparent (0, ii) )] +
                        state_vector [ordUy (quadrant->gparent (1, ii) )]);

      
    }
  }
  
  for (int ii = 0; ii < 4; ++ii){
    const auto& hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    
    vel_rusanov_x[ii] = hpoint>epsilon ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
    vel_rusanov_y[ii] = hpoint>epsilon ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;
    
    const auto dtoptx = hpoint>epsilon ? Dx/vel_rusanov_x[ii] : DELTAT;
    const auto dtopty = hpoint>epsilon ? Dy/vel_rusanov_y[ii] : DELTAT;
    const auto dtopt = dtoptx > dtopty ? dtopty : dtoptx;
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
TG2_scheme::first_step (tmesh::quadrant_iterator quadrant, const double& deltat)
{
  // look at tmesh.h
  const auto & index_quadrant = quadrant->get_global_quad_idx (); // get_forest_quad_idx, get_global_quad_idx
  
  for (int ii = 0; ii < 4; ++ii)
  {
    xn[ii] = quadrant->p (0, ii);
    yn[ii] = quadrant->p (1, ii);
  }
  
  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;
  
  double h_cell_average = 0., Ux_cell_average = 0., Uy_cell_average = 0.;
  double source_h_cell_average = 0., source_Ux_cell_average = 0., source_Uy_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      h_cell_average  += state_vector [ordh (quadrant->gt (ii) )];
      Ux_cell_average += state_vector [ordUx (quadrant->gt (ii) )];
      Uy_cell_average += state_vector [ordUy (quadrant->gt (ii) )];
      
      source_h_cell_average  += h_src_formula  (hdof[ii], Uxdof[ii], Uydof[ii]);
      source_Ux_cell_average += Ux_src_formula (hdof[ii], Uxdof[ii], Uydof[ii]);
      source_Uy_cell_average += Uy_src_formula (hdof[ii], Uxdof[ii], Uydof[ii]);
      
      fluxx_h_node[ii]  = h_flux_formula_x (hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxy_h_node[ii]  = h_flux_formula_y (hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxx_Ux_node[ii] = Ux_flux_formula_x(hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxy_Ux_node[ii] = Ux_flux_formula_y(hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxx_Uy_node[ii] = Uy_flux_formula_x(hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxy_Uy_node[ii] = Uy_flux_formula_y(hdof[ii], Uxdof[ii], Uydof[ii]);
      
    }
    else
    {
      h_cell_average  += .5 * (state_vector [ordh (quadrant->gparent (0, ii) )] +
                        state_vector [ordh (quadrant->gparent (1, ii) )]);
      Ux_cell_average += .5 * (state_vector [ordUx (quadrant->gparent (0, ii) )] +
                        state_vector [ordUx (quadrant->gparent (1, ii) )]);
      Uy_cell_average += .5 * (state_vector [ordUy (quadrant->gparent (0, ii) )] +
                        state_vector [ordUy (quadrant->gparent (1, ii) )]);
      
      source_h_cell_average  += h_src_formula  (hdof[ii], Uxdof[ii], Uydof[ii]);
      source_Ux_cell_average += Ux_src_formula (hdof[ii], Uxdof[ii], Uydof[ii]);//, quadrant);
      source_Uy_cell_average += Uy_src_formula (hdof[ii], Uxdof[ii], Uydof[ii]);//, quadrant);
      
      fluxx_h_node[ii]  = h_flux_formula_x (hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxy_h_node[ii]  = h_flux_formula_y (hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxx_Ux_node[ii] = Ux_flux_formula_x(hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxy_Ux_node[ii] = Ux_flux_formula_y(hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxx_Uy_node[ii] = Uy_flux_formula_x(hdof[ii], Uxdof[ii], Uydof[ii]);
      fluxy_Uy_node[ii] = Uy_flux_formula_y(hdof[ii], Uxdof[ii], Uydof[ii]);
      
    }
  }
  h_cell_average  /= 4.;
  Ux_cell_average /= 4.;
  Uy_cell_average /= 4.;
  
  source_h_cell_average  /= 4.;
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
  
  
  sol_onehalf[index_quadrant]  = h_cell_average - deltat / 2. * div_Fh_cell / area  + deltat / 2. * source_h_cell_average;
  sol_onehalf[index_quadrant+gn_elements] = Ux_cell_average - deltat / 2. * div_FUx_cell / area + deltat / 2. * source_Ux_cell_average;
  sol_onehalf[index_quadrant+2*gn_elements] = Uy_cell_average - deltat / 2. * div_FUy_cell / area + deltat / 2. * source_Uy_cell_average;

}



void
TG2_scheme::second_step (tmesh::quadrant_iterator quadrant, Q1& increment)
{
  
  // look at tmesh.h
  const auto & index_quadrant = quadrant->get_global_quad_idx ();
  
  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }
  
  Dx = xn[1]-xn[0];
  Dy = yn[2]-yn[0];
  area = Dx * Dy;
  
  
  for (int ii = 0; ii < 4; ++ii){
    if (! quadrant->is_hanging (ii)){
      hdof[ii]  = state_vector [ordh (quadrant->gt (ii))];
      Uxdof[ii] = state_vector [ordUx (quadrant->gt (ii))];
      Uydof[ii] = state_vector [ordUy (quadrant->gt (ii))];
      
      isdof_or_hanging[ii] = 1.;
    } else {
      hdof[ii]  = .5 * (state_vector [ordh (quadrant->gparent(0,ii))] +
                        state_vector [ordh (quadrant->gparent(1,ii))]);
      Uxdof[ii] = .5 * (state_vector [ordUx (quadrant->gparent(0,ii))] +
                        state_vector [ordUx (quadrant->gparent(1,ii))]);
      Uydof[ii] = .5 * (state_vector [ordUy (quadrant->gparent(0,ii))] +
                        state_vector [ordUy (quadrant->gparent(1,ii))]);
      
      isdof_or_hanging[ii] = .5;
    }
  }
  
  // weights coefficients for the flux term
  der_coeffs_x = {-Dy/2.*isdof_or_hanging[0], +Dy/2.*isdof_or_hanging[1],
    -Dy/2.*isdof_or_hanging[2], +Dy/2.*isdof_or_hanging[3]};
  
  der_coeffs_y = {-Dx/2.*isdof_or_hanging[0], -Dx/2.*isdof_or_hanging[1],
    +Dx/2.*isdof_or_hanging[2], +Dx/2.*isdof_or_hanging[3]};
  
  double F_star_h_x = 0., F_star_h_y = 0.,
  F_star_Ux_x = 0., F_star_Ux_y = 0.,
  F_star_Uy_x = 0., F_star_Uy_y = 0.;
  
//  for (int ii = 0; ii < 4; ++ii)
//  {
//    F_star_h_x += h_flux_formula_x  (hdof[ii], Uxdof[ii], Uydof[ii]);
//    F_star_h_y += h_flux_formula_y  (hdof[ii], Uxdof[ii], Uydof[ii]);
//
//    F_star_Ux_x += Ux_flux_formula_x (hdof[ii], Uxdof[ii], Uydof[ii]);
//    F_star_Ux_y += Ux_flux_formula_y (hdof[ii], Uxdof[ii], Uydof[ii]);
//
//    F_star_Uy_x += Uy_flux_formula_x (hdof[ii], Uxdof[ii], Uydof[ii]);
//    F_star_Uy_y += Uy_flux_formula_y (hdof[ii], Uxdof[ii], Uydof[ii]);
//  }
//  F_star_h_x  /= 4.;
//  F_star_h_y  /= 4.;
//  F_star_Ux_x /= 4.;
//  F_star_Ux_y /= 4.;
//  F_star_Uy_x /= 4.;
//  F_star_Uy_y /= 4.;
  
  
  F_star_h_x  = h_flux_formula_x(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
  F_star_h_y  = h_flux_formula_y(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);

  F_star_Ux_x = Ux_flux_formula_x(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
  F_star_Ux_y = Ux_flux_formula_y(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);

  F_star_Uy_x = Uy_flux_formula_x(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
  F_star_Uy_y = Uy_flux_formula_y(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
  
  

  
  const double & h_cell  = sol_onehalf[index_quadrant];
  const double & Ux_cell = sol_onehalf[index_quadrant+gn_elements];
  const double & Uy_cell = sol_onehalf[index_quadrant+2*gn_elements];
  
  // compute local extrema
  double h_min  = h_cell,  h_max  = h_cell;
  double Ux_min = Ux_cell, Ux_max = Ux_cell;
  double Uy_min = Uy_cell, Uy_max = Uy_cell;
  for (auto quadrant_nei = quadrant->begin_neighbor_sweep();
       quadrant_nei != quadrant->end_neighbor_sweep (); ++quadrant_nei)
  {
    const auto & index_quadrant_nei = quadrant_nei->get_global_quad_idx ();
    
    const auto & h_current_cell  = local_extrema[index_quadrant_nei];
    const auto & Ux_current_cell = local_extrema[index_quadrant_nei+gn_elements];
    const auto & Uy_current_cell = local_extrema[index_quadrant_nei+2*gn_elements];

    h_min  = std::min(h_current_cell, h_min);
    h_max  = std::max(h_current_cell, h_max);

    Ux_min = std::min(Ux_current_cell, Ux_min);
    Ux_max = std::max(Ux_current_cell, Ux_max);

    Uy_min = std::min(Uy_current_cell, Uy_min);
    Uy_max = std::max(Uy_current_cell, Uy_max);
  }
  
  for (int ii = 0; ii < 4; ++ii){
    const auto& hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    vel_rusanov_cell_x += hpoint>epsilon ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
    vel_rusanov_cell_y += hpoint>epsilon ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;
  }
  vel_rusanov_cell_x /= 4.;
  vel_rusanov_cell_y /= 4.;
  
  
  grad_cell_h  = {.5 * ( (hdof[3]  - hdof[2])  + (hdof[1]  - hdof[0])  ) / Dx, .5 * ( (hdof[2]  - hdof[0])  + (hdof[3] - hdof[1])  ) / Dy};
  grad_cell_Ux = {.5 * ( (Uxdof[3] - Uxdof[2]) + (Uxdof[1] - Uxdof[0]) ) / Dx, .5 * ( (Uxdof[2] - Uxdof[0]) + (Uxdof[3] - Uxdof[1]) ) / Dy};
  grad_cell_Uy = {.5 * ( (Uydof[3] - Uydof[2]) + (Uydof[1] - Uydof[0]) ) / Dx, .5 * ( (Uydof[2] - Uydof[0]) + (Uydof[3] - Uydof[1]) ) / Dy};
  
  // compute flux limiter, grad limiter
  const auto toll = 1e-6;
  double phi_cell_h = 1., phi_cell_Ux = 1., phi_cell_Uy = 1.;
  for (int ii = 0; ii < 4; ++ii){
    const auto h_vertex  = state_vector[ordh   (quadrant->gt (ii))];
    const auto Ux_vertex = state_vector[ordUx  (quadrant->gt (ii))];
    const auto Uy_vertex = state_vector[ordUy  (quadrant->gt (ii))];
    
//    const auto x_node = quadrant->p(0, ii);
//    const auto y_node = quadrant->p(1, ii);
//    const auto x_center = quadrant->centroid(0);
//    const auto y_center = quadrant->centroid(1);
//
//    const auto Dx_v = x_node - x_center;
//    const auto Dy_v = y_node - y_center;
//
//    const auto h_vertex  = h_cell  + grad_cell_h[0] *Dx_v/2 + grad_cell_h[1] *Dy_v/2;
//    const auto Ux_vertex = Ux_cell + grad_cell_Ux[0]*Dx_v/2 + grad_cell_Ux[1]*Dy_v/2;
//    const auto Uy_vertex = Uy_cell + grad_cell_Uy[0]*Dx_v/2 + grad_cell_Uy[1]*Dy_v/2;
    
    flux_limiter(h_min,  h_max,  h_vertex,  h_cell,  toll, phi_cell_h);
    flux_limiter(Ux_min, Ux_max, Ux_vertex, Ux_cell, toll, phi_cell_Ux);
    flux_limiter(Uy_min, Uy_max, Uy_vertex, Uy_cell, toll, phi_cell_Uy);
    
  }
  phi_cell_h  = 1. - phi_cell_h;
  phi_cell_Ux = 1. - phi_cell_Ux;
  phi_cell_Uy = 1. - phi_cell_Uy;
  
  
  F_star_h_x  -= grad_cell_h[0]*Dx  * phi_cell_h  * vel_rusanov_cell_y;
  F_star_h_y  -= grad_cell_h[1]*Dy  * phi_cell_h  * vel_rusanov_cell_x;
  
  F_star_Ux_x -= grad_cell_Ux[0]*Dx * phi_cell_Ux * vel_rusanov_cell_y;
  F_star_Ux_y -= grad_cell_Ux[1]*Dy * phi_cell_Ux * vel_rusanov_cell_x;
  
  F_star_Uy_x -= grad_cell_Uy[0]*Dx * phi_cell_Uy * vel_rusanov_cell_y;
  F_star_Uy_y -= grad_cell_Uy[1]*Dy * phi_cell_Uy * vel_rusanov_cell_x;
  
  
  
  
  
  
  
  
  // aggiungere termine sorgente!!
  for (int ii = 0; ii < 4; ++ii){
    
    const auto h_  = der_coeffs_x[ii] * F_star_h_x  + der_coeffs_y[ii] * F_star_h_y;
    const auto Ux_ = der_coeffs_x[ii] * F_star_Ux_x + der_coeffs_y[ii] * F_star_Ux_y;
    const auto Uy_ = der_coeffs_x[ii] * F_star_Uy_x + der_coeffs_y[ii] * F_star_Uy_y;
    
    if (! quadrant->is_hanging (ii)){
      increment [ordh  (quadrant->gt (ii))] += h_;
      increment [ordUx (quadrant->gt (ii))] += Ux_;
      increment [ordUy (quadrant->gt (ii))] += Uy_;
      
      
      auto boundary_idxx = quadrant->ex (ii);
      auto boundary_idxy = quadrant->ey (ii);
      
      
      
      if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY) {
        double F_star_h_x_b = 0., F_star_h_y_b = 0.,
        F_star_Ux_x_b = 0., F_star_Ux_y_b = 0., F_star_Uy_x_b = 0., F_star_Uy_y_b = 0.;
        
        F_star_h_x_b = h_flux_formula_x(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
        F_star_h_y_b = h_flux_formula_y(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);

        F_star_Ux_x_b = Ux_flux_formula_x(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
        F_star_Ux_y_b = Ux_flux_formula_y(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);

        F_star_Uy_x_b = Uy_flux_formula_x(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
        F_star_Uy_y_b = Uy_flux_formula_y(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]); 
    
        
        
        auto h_b  = -der_coeffs_x[ii] * F_star_h_x_b  + der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = -der_coeffs_x[ii] * F_star_Ux_x_b + der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = -der_coeffs_x[ii] * F_star_Uy_x_b + der_coeffs_y[ii] * F_star_Uy_y_b;
        
//        auto h_b  = -der_coeffs_x[ii] * F_star_h_x  + der_coeffs_y[ii] * F_star_h_y;
//        auto Ux_b = -der_coeffs_x[ii] * F_star_Ux_x + der_coeffs_y[ii] * F_star_Ux_y;
//        auto Uy_b = -der_coeffs_x[ii] * F_star_Uy_x + der_coeffs_y[ii] * F_star_Uy_y;
        
        
        increment [ordh  (quadrant->gt (ii))] += h_b;
        increment [ordUx (quadrant->gt (ii))] += Ux_b;
        increment [ordUy (quadrant->gt (ii))] += Uy_b;
        
        
        for (int jj = 0; jj < 4; ++jj){
          if ( quadrant->is_hanging (jj) && (der_coeffs_y[jj]*der_coeffs_y[ii] > 0.) )
          {
            F_star_h_y_b = h_flux_formula_y(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
            
            F_star_Ux_y_b = Ux_flux_formula_y(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
            
            F_star_Uy_y_b = Uy_flux_formula_y(sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant+gn_elements], sol_onehalf[index_quadrant+2*gn_elements]);
            
            h_b  = 2. * der_coeffs_y[jj] * F_star_h_y_b;
            Ux_b = 2. * der_coeffs_y[jj] * F_star_Ux_y_b;
            Uy_b = 2. * der_coeffs_y[jj] * F_star_Uy_y_b;
            
//            h_b  = 2. * der_coeffs_y[jj] * F_star_h_y;
//            Ux_b = 2. * der_coeffs_y[jj] * F_star_Ux_y;
//            Uy_b = 2. * der_coeffs_y[jj] * F_star_Uy_y;
            
            increment [ordh  (quadrant->gt (ii))] += h_b;
            increment [ordUx (quadrant->gt (ii))] += Ux_b;
            increment [ordUy (quadrant->gt (ii))] += Uy_b;
          }
          
        }
        
        
        
        
      }
      
      if (boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY) {
        double F_star_h_x_b = 0., F_star_h_y_b = 0.,
        F_star_Ux_x_b = 0., F_star_Ux_y_b = 0., F_star_Uy_x_b = 0., F_star_Uy_y_b = 0.;
        
        F_star_h_x_b = h_flux_formula_x(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], -sol_onehalf[index_quadrant+2*gn_elements]);
        F_star_h_y_b = h_flux_formula_y(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], -sol_onehalf[index_quadrant+2*gn_elements]);

        F_star_Ux_x_b = Ux_flux_formula_x(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], -sol_onehalf[index_quadrant+2*gn_elements]);
        F_star_Ux_y_b = Ux_flux_formula_y(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], -sol_onehalf[index_quadrant+2*gn_elements]);

        F_star_Uy_x_b = Uy_flux_formula_x(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], -sol_onehalf[index_quadrant+2*gn_elements]);
        F_star_Uy_y_b = Uy_flux_formula_y(sol_onehalf[index_quadrant], sol_onehalf[index_quadrant+gn_elements], -sol_onehalf[index_quadrant+2*gn_elements]);
        
        
        auto h_b  = der_coeffs_x[ii] * F_star_h_x_b  - der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = der_coeffs_x[ii] * F_star_Ux_x_b - der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = der_coeffs_x[ii] * F_star_Uy_x_b - der_coeffs_y[ii] * F_star_Uy_y_b;
        
//        auto h_b  = der_coeffs_x[ii] * F_star_h_x  - der_coeffs_y[ii] * F_star_h_y;
//        auto Ux_b = der_coeffs_x[ii] * F_star_Ux_x - der_coeffs_y[ii] * F_star_Ux_y;
//        auto Uy_b = der_coeffs_x[ii] * F_star_Uy_x - der_coeffs_y[ii] * F_star_Uy_y;
        
        
        increment [ordh  (quadrant->gt (ii))] += h_b;
        increment [ordUx (quadrant->gt (ii))] += Ux_b;
        increment [ordUy (quadrant->gt (ii))] += Uy_b;
        
        
        for (int jj = 0; jj < 4; ++jj){
          if ( quadrant->is_hanging (jj) && (der_coeffs_x[jj]*der_coeffs_x[ii] > 0.) )
          {
            F_star_h_x_b = h_flux_formula_x (sol_onehalf[index_quadrant], sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
            
            F_star_Ux_x_b = Ux_flux_formula_x (sol_onehalf[index_quadrant], sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
            
            F_star_Uy_x_b = Uy_flux_formula_x (sol_onehalf[index_quadrant], sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
            
            h_b  = 2. * der_coeffs_x[jj] * F_star_h_x_b;
            Ux_b = 2. * der_coeffs_x[jj] * F_star_Ux_x_b;
            Uy_b = 2. * der_coeffs_x[jj] * F_star_Uy_x_b;
            
//            h_b  = 2. * der_coeffs_x[jj] * F_star_h_x;
//            Ux_b = 2. * der_coeffs_x[jj] * F_star_Ux_x;
//            Uy_b = 2. * der_coeffs_x[jj] * F_star_Uy_x;
            
            increment [ordh  (quadrant->gt (ii))] += h_b;
            increment [ordUx (quadrant->gt (ii))] += Ux_b;
            increment [ordUy (quadrant->gt (ii))] += Uy_b;
          }
          
        }
        
        
      }
      
      // corner points!
      if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY && boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY)
      {
        double F_star_h_x_b = 0., F_star_h_y_b = 0.,
        F_star_Ux_x_b = 0., F_star_Ux_y_b = 0., F_star_Uy_x_b = 0., F_star_Uy_y_b = 0.;
        
        F_star_h_x_b = h_flux_formula_x (sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
        F_star_h_y_b = h_flux_formula_y (sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
        
        F_star_Ux_x_b = Ux_flux_formula_x (sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
        F_star_Ux_y_b = Ux_flux_formula_y (sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
        
        F_star_Uy_x_b = Uy_flux_formula_x (sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
        F_star_Uy_y_b = Uy_flux_formula_y (sol_onehalf[index_quadrant], -sol_onehalf[index_quadrant + gn_elements], -sol_onehalf[index_quadrant + 2 * gn_elements]);
        
        
        auto h_b  = -der_coeffs_x[ii] * F_star_h_x_b  - der_coeffs_y[ii] * F_star_h_y_b;
        auto Ux_b = -der_coeffs_x[ii] * F_star_Ux_x_b - der_coeffs_y[ii] * F_star_Ux_y_b;
        auto Uy_b = -der_coeffs_x[ii] * F_star_Uy_x_b - der_coeffs_y[ii] * F_star_Uy_y_b;
        
//        auto h_b  = -der_coeffs_x[ii] * F_star_h_x  - der_coeffs_y[ii] * F_star_h_y;
//        auto Ux_b = -der_coeffs_x[ii] * F_star_Ux_x - der_coeffs_y[ii] * F_star_Ux_y;
//        auto Uy_b = -der_coeffs_x[ii] * F_star_Uy_x - der_coeffs_y[ii] * F_star_Uy_y;
        
        
        increment [ordh  (quadrant->gt (ii))] += h_b;
        increment [ordUx (quadrant->gt (ii))] += Ux_b;
        increment [ordUy (quadrant->gt (ii))] += Uy_b;
      
        
      }
      
      
    } else {
      
      increment [ordh  (quadrant->gparent(0,ii))] += h_;
      increment [ordh  (quadrant->gparent(1,ii))] += h_;
      
      increment [ordUx (quadrant->gparent(0,ii))] += Ux_;
      increment [ordUx (quadrant->gparent(1,ii))] += Ux_;
      
      increment [ordUy (quadrant->gparent(0,ii))] += Uy_;
      increment [ordUy (quadrant->gparent(1,ii))] += Uy_;
      
      
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
{ return h > epsilon ?  Ux*Ux/h + grav*h*h/2. : 0.; }

double
TG2_scheme::Ux_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ? Uy*Ux/h : 0.; }

double
TG2_scheme::Uy_flux_formula_x (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ? Uy*Ux/h : 0.; }

double
TG2_scheme::Uy_flux_formula_y (const double& h, const double& Ux, const double& Uy)
{ return h > epsilon ? Uy*Uy/h + grav*h*h/2. : 0.; }


// source terms
double
TG2_scheme::h_src_formula (const double& h, const double& Ux, const double& Uy)
{ return (0.); }

double
TG2_scheme::Ux_src_formula (const double& h, const double& Ux, const double& Uy)
{
  const double dZdx = 0.;
  return (-grav*h*dZdx);
}

double
TG2_scheme::Uy_src_formula (const double& h, const double& Ux, const double& Uy)
{
  const double dZdy = 0.;
  return (-grav*h*dZdy);
}



