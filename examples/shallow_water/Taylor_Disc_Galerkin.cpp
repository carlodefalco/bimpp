#include "Taylor_Disc_Galerkin.h"
#include <algorithm>
#include <cassert>

TG2_scheme::TG2_scheme(const distributed_vector& dg_coefficients,
                       const distributed_vector& dg_coefficients_old,
                       const distributed_vector& dg_coefficients_oldold,
                       distributed_vector& incr_dg_coefficients,
                       Q0& sol_onehalf,
                       const ordering& oh,
                       const ordering& oUx,
                       const ordering& oUy,
                       const ordering& oh_nodal,
                       const ordering& oUx_nodal,
                       const ordering& oUy_nodal,
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
: dg_coefficients(dg_coefficients), dg_coefficients_old(dg_coefficients_old), dg_coefficients_oldold(dg_coefficients_oldold), incr_dg_coefficients(incr_dg_coefficients), sol_onehalf(sol_onehalf),
  ordh(oh), ordUx(oUx), ordUy(oUy), ordh_nodal(oh_nodal), ordUx_nodal(oUx_nodal), ordUy_nodal(oUy_nodal), Z(Z), slope_x(slope_x), slope_y(slope_y), DELTAT(DELTAT), epsilon(h_min), is_non_reflBC(is_non_reflBC), is_bed_friction(is_bed_friction), is_stress_tensor(is_stress_tensor), grav(grav),
  density(density), turbulence_coeff(turbulence_coeff), surface_pressure(surface_pressure), bed_friction_angle_rad(bed_friction_angle_rad), fluid_viscosity(fluid_viscosity), yield_shear_stress(yield_shear_stress)
{ }

 

void
TG2_scheme::compute_dt (tmesh::quadrant_iterator quadrant)
{

  Dx = quadrant->p (0, 1) - quadrant->p (0, 0);
  Dy = quadrant->p (1, 2) - quadrant->p (1, 0);
  

  for (int ii = 0; ii < 4; ++ii)
  {
    if (! quadrant->is_hanging (ii) )
    {
      hdof [ii] = dg_coefficients[ordh  (quadrant->get_global_quad_idx ()) + ii];
      Uxdof[ii] = dg_coefficients[ordUx (quadrant->get_global_quad_idx ()) + ii];
      Uydof[ii] = dg_coefficients[ordUy (quadrant->get_global_quad_idx ()) + ii]; 
    }
    else
    {
      hdof [ii] = .5 * dg_coefficients[ordh  (quadrant->get_global_quad_idx ()) + ii];
      Uxdof[ii] = .5 * dg_coefficients[ordUx (quadrant->get_global_quad_idx ()) + ii];
      Uydof[ii] = .5 * dg_coefficients[ordUy (quadrant->get_global_quad_idx ()) + ii];
    }
  }

  
  for (int ii = 0; ii < 4; ++ii){
    const auto & hpoint = hdof[ii];
    const auto celerity = std::sqrt(grav*hpoint);
    
    const auto vel_rusanov_cell_x = hpoint>epsilon ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
    const auto vel_rusanov_cell_y = hpoint>epsilon ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;

    //if (Dx/vel_rusanov_cell_x < 1e-3)
    //std::cout << vel_rusanov_cell_x << " " << celerity << " " << Dx/vel_rusanov_cell_x << " " << Uxdof[ii] << std::endl;

    const auto dtoptx = hpoint>epsilon ? Dx/vel_rusanov_cell_x : DELTAT;
    const auto dtopty = hpoint>epsilon ? Dy/vel_rusanov_cell_y : DELTAT;

    const auto dt_c = std::min(dtoptx, dtopty);
    dt = dt_c < dt ? dt_c : dt;

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
      hdof      = dg_coefficients       [ordh (quadrant->get_global_quad_idx ()) + ii]; 
      hdof_old  = dg_coefficients_old   [ordh (quadrant->get_global_quad_idx ()) + ii]; 
      hdof_oldd = dg_coefficients_oldold[ordh (quadrant->get_global_quad_idx ()) + ii];
    }
    else
    {
      hdof      = .5 * dg_coefficients       [ordh (quadrant->get_global_quad_idx ()) + ii];
      hdof_old  = .5 * dg_coefficients_old   [ordh (quadrant->get_global_quad_idx ()) + ii];
      hdof_oldd = .5 * dg_coefficients_oldold[ordh (quadrant->get_global_quad_idx ()) + ii];
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
  
  Dx = quadrant->p (0, 1) - quadrant->p (0, 0);
  Dy = quadrant->p (1, 2) - quadrant->p (1, 0);
  area = Dx * Dy;


  
  double h_cell_average = 0., Ux_cell_average = 0., Uy_cell_average = 0., source_Ux_cell_average = 0., source_Uy_cell_average = 0.;
  for (int ii = 0; ii < 4; ++ii)
  {
    double hdof_c, Uxdof_c, Uydof_c;

    if (! quadrant->is_hanging (ii) )
    {
      hdof_c    = dg_coefficients[ordh  (quadrant->get_global_quad_idx ()) + ii];
      Uxdof_c   = dg_coefficients[ordUx (quadrant->get_global_quad_idx ()) + ii];
      Uydof_c   = dg_coefficients[ordUy (quadrant->get_global_quad_idx ()) + ii];
    }
    else
    {
      hdof_c    = .5 * dg_coefficients[ordh  (quadrant->get_global_quad_idx ()) + ii];
      Uxdof_c   = .5 * dg_coefficients[ordUx (quadrant->get_global_quad_idx ()) + ii];
      Uydof_c   = .5 * dg_coefficients[ordUy (quadrant->get_global_quad_idx ()) + ii];
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

  //std::cout << h_cell_average  - dt / 2. * div_Fh_cell  / area << std::endl;
  
  sol_onehalf[ordh_nodal    (index_quadrant)] = h_cell_average  - dt/2. * div_Fh_cell  / area;
  sol_onehalf[ordUx_nodal   (index_quadrant)] = Ux_cell_average - dt/2. * div_FUx_cell / area + dt/2. * source_Ux_cell_average;
  sol_onehalf[ordUy_nodal   (index_quadrant)] = Uy_cell_average - dt/2. * div_FUy_cell / area + dt/2. * source_Uy_cell_average;

}


void
TG2_scheme::second_step (tmesh::quadrant_iterator quadrant)
{

  // look at tmesh.h
  const auto & index_quadrant = quadrant->get_forest_quad_idx (); 

  std::array<int,4> bimpp_to_rev_ord = {0, 1, 3, 2};

  for (int ii = 0; ii < 4; ++ii) {
    xn[ii] = quadrant->p(0, ii);
    yn[ii] = quadrant->p(1, ii);
  }

  Dx = xn[1] - xn[0];
  Dy = yn[2] - yn[0];
  area = Dx * Dy;

  for (int ii = 0; ii < 4; ++ii){
      if (! quadrant->is_hanging (ii)){
        isdof_or_hanging[ii] = 1.;
      } else {
        isdof_or_hanging[ii] = .5;
      }
  }


  // compute boundary fluxes --> bisogna mettere anche i flussi degli sforzi!!!!!!!!!!!!!!!!
  bool is_boundary_edge = true;

  const auto & h_cell  = sol_onehalf[ordh_nodal    (index_quadrant)];
  const auto & Ux_cell = sol_onehalf[ordUx_nodal   (index_quadrant)];
  const auto & Uy_cell = sol_onehalf[ordUy_nodal   (index_quadrant)];


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

      const auto & h_cell_nei  = sol_onehalf[ordh_nodal   (index_quadrant_nei)];
      const auto & Ux_cell_nei = sol_onehalf[ordUx_nodal  (index_quadrant_nei)];
      const auto & Uy_cell_nei = sol_onehalf[ordUy_nodal  (index_quadrant_nei)];


      for (int jEdge = 0; jEdge < 4; ++jEdge) { // cycle neigh edges 

        const auto j_1 = bimpp_to_rev_ord[jEdge];
        const auto j_2 = bimpp_to_rev_ord[(jEdge+1)%4];

        const auto edge_length_nei = std::sqrt(std::pow((Xn[j_1]-Xn[j_2]),2.) + std::pow((Yn[j_1]-Yn[j_2]),2.));
        const std::array<double,2> outward_normal_edge_nei = {(-Yn[j_1]+Yn[j_2])/edge_length_nei, ( Xn[j_1]-Xn[j_2])/edge_length_nei};  
        //const std::array<double,2> array_edge = {(Xn[j_2]-Xn[j_1])/edge_length_nei, (Yn[j_2]-Yn[j_1])/edge_length_nei};
        const bool check_orthogonality = std::inner_product(outward_normal_edge_nei.begin(), outward_normal_edge_nei.end(), outward_normal_edge.begin(), 0.) == -1;

        //std::cout << check_orthogonality << std::endl;
        
        // std::cout << index_quadrant << " " << index_quadrant_nei << " " << ((((xn[i_1] == Xn[j_1] && yn[i_1] == Yn[j_1]) || 
        //        (xn[i_2] == Xn[j_1] && yn[i_2] == Yn[j_1]))||
        //       ((xn[i_1] == Xn[j_2] && yn[i_1] == Yn[j_2]) || 
        //        (xn[i_2] == Xn[j_2] && yn[i_2] == Yn[j_2]))) && check_orthogonality && index_quadrant!=index_quadrant_nei) << " " << (((xn[i_1] == Xn[j_1] && yn[i_1] == Yn[j_1]) || 
        //        (xn[i_2] == Xn[j_1] && yn[i_2] == Yn[j_1]))||
        //       ((xn[i_1] == Xn[j_2] && yn[i_1] == Yn[j_2]) || 
        //        (xn[i_2] == Xn[j_2] && yn[i_2] == Yn[j_2]))) << " " << check_orthogonality << " " << iEdge << " " << is_boundary_edge << std::endl;

        // std::cout << quadrant->centroid(0) << " " << quadrant->centroid(1) << " " << quadrant_nei->centroid(0) << " " << quadrant_nei->centroid(1) << " " << ( (((xn[i_1] == Xn[j_1] && yn[i_1] == Yn[j_1]) || 
        //        (xn[i_2] == Xn[j_1] && yn[i_2] == Yn[j_1]))||
        //       ((xn[i_1] == Xn[j_2] && yn[i_1] == Yn[j_2]) || 
        //        (xn[i_2] == Xn[j_2] && yn[i_2] == Yn[j_2]))) && check_orthogonality && index_quadrant!=index_quadrant_nei ) << " " << iEdge << " " << index_quadrant << " " << index_quadrant_nei << std::endl;// << " " << outward_normal_edge_nei[0] << " " << outward_normal_edge_nei[1] << " " << outward_normal_edge[0] << " " << outward_normal_edge[1] << " " << std::inner_product(outward_normal_edge_nei.begin(), outward_normal_edge_nei.end(), outward_normal_edge.begin(), 0.) << std::endl;

        if ( (((xn[i_1] == Xn[j_1] && yn[i_1] == Yn[j_1]) || 
               (xn[i_2] == Xn[j_1] && yn[i_2] == Yn[j_1]))||
              ((xn[i_1] == Xn[j_2] && yn[i_1] == Yn[j_2]) || 
               (xn[i_2] == Xn[j_2] && yn[i_2] == Yn[j_2]))) && check_orthogonality && index_quadrant!=index_quadrant_nei )
        {

          is_boundary_edge = false;


          const auto edge_boundary = std::min(edge_length, edge_length_nei);

          const auto speed     = h_cell    >epsilon ? std::abs((Ux_cell    /h_cell    )*outward_normal_edge[0] + (Uy_cell    /h_cell    )*outward_normal_edge[1]) + std::sqrt(grav*h_cell    ) : 0.;
          const auto speed_nei = h_cell_nei>epsilon ? std::abs((Ux_cell_nei/h_cell_nei)*outward_normal_edge[0] + (Uy_cell_nei/h_cell_nei)*outward_normal_edge[1]) + std::sqrt(grav*h_cell_nei) : 0.;

          const auto smax = std::max(speed, speed_nei);

          // interface Rusanov fluxes
          const auto flux_int_h  = .5*((h_flux_formula_x (h_cell, Ux_cell, Uy_cell)+h_flux_formula_x (h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (h_flux_formula_y (h_cell, Ux_cell, Uy_cell)+h_flux_formula_y (h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(h_cell_nei -h_cell );
          const auto flux_int_Ux = .5*((Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(Ux_cell_nei-Ux_cell);
          const auto flux_int_Uy = .5*((Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uy_cell_nei-Uy_cell);

          //std::cout << edge_boundary << " " << flux_int_h << " " << flux_int_Ux << " " << flux_int_Uy << std::endl;// " " << outward_normal_edge[0] << " " << outward_normal_edge[1] << " " << smax << " " << h_cell_nei << std::endl;

          //std::cout << flux_int_h << " " << incr_dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + i_1]; //quadrant->get_global_quad_idx () << " " << ordh   (quadrant->get_global_quad_idx ()) + i_1 << " " << ordh   (quadrant->get_global_quad_idx ()) + i_2 << std::endl;

          // occhio al segno qui!!!!!!!!!!!!!!!!
          incr_dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + i_1] += -edge_boundary*flux_int_h *.5*isdof_or_hanging[i_1];
          incr_dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + i_2] += -edge_boundary*flux_int_h *.5*isdof_or_hanging[i_2];

          incr_dg_coefficients[ordUx  (quadrant->get_global_quad_idx ()) + i_1] += -edge_boundary*flux_int_Ux*.5*isdof_or_hanging[i_1];
          incr_dg_coefficients[ordUx  (quadrant->get_global_quad_idx ()) + i_2] += -edge_boundary*flux_int_Ux*.5*isdof_or_hanging[i_2];

          incr_dg_coefficients[ordUy  (quadrant->get_global_quad_idx ()) + i_1] += -edge_boundary*flux_int_Uy*.5*isdof_or_hanging[i_1];
          incr_dg_coefficients[ordUy  (quadrant->get_global_quad_idx ()) + i_2] += -edge_boundary*flux_int_Uy*.5*isdof_or_hanging[i_2];
          
          //std::cout << " " << incr_dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + i_1] << std::endl;
          
          break; // this just goes outside the jEdge cycle 
        }
      }

    }

    //if (is_boundary_edge)
    //std::cout << is_boundary_edge << std::endl;

    if (is_boundary_edge) // set boundary conditions
    { 
            
      auto h_cell_nei  = h_cell;
      auto Ux_cell_nei = Ux_cell;
      auto Uy_cell_nei = Uy_cell;

      //std::cout << (!is_non_reflBC) << std::endl;

      Ux_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell + outward_normal_edge[1]*Uy_cell)*outward_normal_edge[0];
      Uy_cell_nei -= (!is_non_reflBC)*2.*(outward_normal_edge[0]*Ux_cell + outward_normal_edge[1]*Uy_cell)*outward_normal_edge[1];

      const auto speed     = h_cell    >epsilon ? std::abs((Ux_cell    /h_cell    )*outward_normal_edge[0] + (Uy_cell    /h_cell    )*outward_normal_edge[1]) + std::sqrt(grav*h_cell    ) : 0.;
      const auto speed_nei = h_cell_nei>epsilon ? std::abs((Ux_cell_nei/h_cell_nei)*outward_normal_edge[0] + (Uy_cell_nei/h_cell_nei)*outward_normal_edge[1]) + std::sqrt(grav*h_cell_nei) : 0.;

      const auto smax = std::max(speed, speed_nei); 

      const auto flux_int_h  = .5*((h_flux_formula_x (h_cell, Ux_cell, Uy_cell)+h_flux_formula_x (h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (h_flux_formula_y (h_cell, Ux_cell, Uy_cell)+h_flux_formula_y (h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(h_cell_nei -h_cell );
      const auto flux_int_Ux = .5*((Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Ux_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(Ux_cell_nei-Ux_cell);
      const auto flux_int_Uy = .5*((Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_x(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[0] + (Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell)+Uy_flux_formula_y(h_cell_nei, Ux_cell_nei, Uy_cell_nei))*outward_normal_edge[1]) - .5*smax*(Uy_cell_nei-Uy_cell);

      //std::cout << flux_int_h << " " << flux_int_Ux << " " << flux_int_Uy << " " << h_cell << " " << grav*h_cell*h_cell/2. << std::endl;
      
      // occhio al segno qui!!!!!!!!!!!!!!!!
      incr_dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + i_1] += -edge_length*flux_int_h *.5*isdof_or_hanging[i_1];
      incr_dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + i_2] += -edge_length*flux_int_h *.5*isdof_or_hanging[i_2];

      incr_dg_coefficients[ordUx  (quadrant->get_global_quad_idx ()) + i_1] += -edge_length*flux_int_Ux*.5*isdof_or_hanging[i_1];
      incr_dg_coefficients[ordUx  (quadrant->get_global_quad_idx ()) + i_2] += -edge_length*flux_int_Ux*.5*isdof_or_hanging[i_2];

      incr_dg_coefficients[ordUy  (quadrant->get_global_quad_idx ()) + i_1] += -edge_length*flux_int_Uy*.5*isdof_or_hanging[i_1];
      incr_dg_coefficients[ordUy  (quadrant->get_global_quad_idx ()) + i_2] += -edge_length*flux_int_Uy*.5*isdof_or_hanging[i_2];
      

    }

  }




  // weight coefficients for the flux term
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
  

  const auto diff_term_h_x  = grad_cell_h [0] * vel_rusanov_cell_y;
  const auto diff_term_h_y  = grad_cell_h [1] * vel_rusanov_cell_x;

  const auto diff_term_Ux_x = grad_cell_Ux[0] * vel_rusanov_cell_y;
  const auto diff_term_Ux_y = grad_cell_Ux[1] * vel_rusanov_cell_x;

  const auto diff_term_Uy_x = grad_cell_Uy[0] * vel_rusanov_cell_y;
  const auto diff_term_Uy_y = grad_cell_Uy[1] * vel_rusanov_cell_x;
  


  // const auto F_star_h_x  = h_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_h_x;
  // const auto F_star_h_y  = h_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_h_y;

  // const auto F_star_Ux_x = Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_x;
  // const auto F_star_Ux_y = Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Ux_y;

  // const auto F_star_Uy_x = Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_x;
  // const auto F_star_Uy_y = Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell) - diff_term_Uy_y;



  const auto F_star_h_x  = h_flux_formula_x(h_cell, Ux_cell, Uy_cell);
  const auto F_star_h_y  = h_flux_formula_y(h_cell, Ux_cell, Uy_cell);

  const auto F_star_Ux_x = Ux_flux_formula_x(h_cell, Ux_cell, Uy_cell);
  const auto F_star_Ux_y = Ux_flux_formula_y(h_cell, Ux_cell, Uy_cell);

  const auto F_star_Uy_x = Uy_flux_formula_x(h_cell, Ux_cell, Uy_cell);
  const auto F_star_Uy_y = Uy_flux_formula_y(h_cell, Ux_cell, Uy_cell);


  // element integrals
  for (int ii = 0; ii < 4; ++ii){
    
    sigma_stress = is_stress_tensor ? compute_nodal_stress (hdof[ii], Uxdof[ii], Uydof[ii], grad_cell_ux, grad_cell_uy) : sigma_stress;

    const auto D_Ux_x = Ux_stress_formula_x(hdof[ii], Uxdof[ii], Uydof[ii]);
    const auto D_Ux_y = Ux_stress_formula_y(hdof[ii], Uxdof[ii], Uydof[ii]);

    const auto D_Uy_x = Uy_stress_formula_x(hdof[ii], Uxdof[ii], Uydof[ii]);
    const auto D_Uy_y = Uy_stress_formula_y(hdof[ii], Uxdof[ii], Uydof[ii]);

    //std::cout << D_Ux_x << " " << D_Ux_y << " " << D_Uy_x << " " << D_Uy_y << std::endl;


    const auto h_  = der_coeffs_x[ii]*F_star_h_x +der_coeffs_y[ii]*F_star_h_y;
    const auto Ux_ = der_coeffs_x[ii]*F_star_Ux_x+der_coeffs_y[ii]*F_star_Ux_y;// + der_coeffs_x[ii]*(1./3.)*D_Ux_x+der_coeffs_y[ii]*(1./3.)*D_Ux_y + .25*area*isdof_or_hanging[ii]*Ux_src_formula(h_cell, Ux_cell, Uy_cell, slope_x[index_quadrant]);
    const auto Uy_ = der_coeffs_x[ii]*F_star_Uy_x+der_coeffs_y[ii]*F_star_Uy_y;// + der_coeffs_x[ii]*(1./3.)*D_Uy_x+der_coeffs_y[ii]*(1./3.)*D_Uy_y + .25*area*isdof_or_hanging[ii]*Uy_src_formula(h_cell, Ux_cell, Uy_cell, slope_y[index_quadrant]);

    
    incr_dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + ii] += h_;
    incr_dg_coefficients[ordUx  (quadrant->get_global_quad_idx ()) + ii] += Ux_;
    incr_dg_coefficients[ordUy  (quadrant->get_global_quad_idx ()) + ii] += Uy_;
    

    //std::cout << "here!!! " << h_ << " " << Ux_ << " " << Uy_ << " " << incr_dg_coefficients[ordh   (quadrant->get_global_quad_idx ()) + ii] << std::endl;

  }



}




void
TG2_scheme::minmod_modified_tvb(const double& a, const double& b, const double& c, const double& Dx, const double& K)
{

}


double
TG2_scheme::minmod(const double& a, const double& b, const double& c)
{
  const double m = (1./3.)*(signum(a)+signum(b)+signum(c)) * std::min(std::min(std::abs(a),std::abs(b)), std::abs(c));

  return m;
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

//std::cout << viscos << std::endl;//def_grad[0] << " " << def_grad[1] << " " << def_grad[2] << std::endl;
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


double
signum (const double& x)
{ return ((x > 0) ? 1.0 : (x < 0) ? -1.0 : 0.0); }

