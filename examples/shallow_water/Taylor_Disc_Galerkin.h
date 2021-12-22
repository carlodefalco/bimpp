#ifndef TAYLOR_GALERKIN_H
#define TAYLOR_GALERKIN_H

#include <numeric>
#include <bim_distributed_vector.h>
#include <tmesh.h>
#include <quad_operators.h>

 


class TG2_scheme
{
  using Q1  = q1_vec<distributed_vector>;
  using Q0  = std::vector<double>;
  
public:
  
  TG2_scheme(const distributed_vector& dg_coefficients,
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
             const double& yield_shear_stress);
  
  TG2_scheme() = delete;
  
  ~TG2_scheme() = default;
  
  
  void
  compute_dt (tmesh::quadrant_iterator quadrant);
  
  void
  compute_dt_adaptive (tmesh::quadrant_iterator quadrant);
  
  void
  first_step (tmesh::quadrant_iterator quadrant);
  
  void
  second_step (tmesh::quadrant_iterator quadrant);
  
  void
  minmod_modified_tvb(const double& a, const double& b, const double& c, const double& Dx, const double& K);
  
  double
  minmod(const double& a, const double& b, const double& c);
  
  
  void
  set_dt (const double dt_);
  
  void
  set_times(const double& time, const double& time_old, const double& time_oldd);
  
  double
  get_dt ();
  
  
  double dt;
  
  double Dx, Dy, area;
  
  
    ///  quadrant vertex (dofs) coordinates
    ///  The assumed numbering for quadrant nodes is
    ///  the following :
    ///    ^
    ///   yI
    ///   2------------------3
    ///   |                  |
    ///   |                  |
    ///   |                  |
    ///   |                  |
    ///   0------------------1 -->x
  
  std::array<double, 4> xn = {0, 0, 0, 0};
  std::array<double, 4> yn = {0, 0, 0, 0};
  
  // local dofs for state vector components
  std::array<double, 4> hdof    = {0, 0, 0, 0};
  std::array<double, 4> Uxdof   = {0, 0, 0, 0};
  std::array<double, 4> Uydof   = {0, 0, 0, 0};
  std::array<double, 4> Z_node  = {0, 0, 0, 0};
  
  // std::array<double, 4> source_h_node  = {0, 0, 0, 0};
  // std::array<double, 4> source_Ux_node = {0, 0, 0, 0};
  // std::array<double, 4> source_Uy_node = {0, 0, 0, 0};
  
  std::array<double, 4> fluxx_h_node    = {0, 0, 0, 0}, fluxy_h_node    = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Ux_node   = {0, 0, 0, 0}, fluxy_Ux_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uy_node   = {0, 0, 0, 0}, fluxy_Uy_node   = {0, 0, 0, 0};
  
  
  std::array<double, 3> sigma_stress = {0., 0., 0.};
  
  
  // flux functions
  double
  h_flux_formula_x (const double& h, const double& Ux, const double& Uy);
  
  double
  h_flux_formula_y (const double& h, const double& Ux, const double& Uy);
  
  double
  Ux_flux_formula_x (const double& h, const double& Ux, const double& Uy);
  
  double
  Ux_flux_formula_y (const double& h, const double& Ux, const double& Uy);
  
  double
  Uy_flux_formula_x (const double& h, const double& Ux, const double& Uy);
  
  double
  Uy_flux_formula_y (const double& h, const double& Ux, const double& Uy);


  // stress functions
  double
  Ux_stress_formula_x (const double& h, const double& Ux, const double& Uy);
  
  double
  Ux_stress_formula_y (const double& h, const double& Ux, const double& Uy);
  
  double
  Uy_stress_formula_x (const double& h, const double& Ux, const double& Uy);
  
  double
  Uy_stress_formula_y (const double& h, const double& Ux, const double& Uy);


  std::array<double,3>
  compute_nodal_stress (const double& h, const double& Ux, const double& Uy, const std::array<double,2>& grad_cell_ux, const std::array<double,2>& grad_cell_uy);

  std::array<double,6>
  compute_nodal_def_grad (const double& h, const double& Ux, const double& Uy, const std::array<double,2>& grad_cell_ux, const std::array<double,2>& grad_cell_uy);


  
  // source terms
  double
  h_src_formula (const double& h, const double& Ux, const double& Uy);
  
  double
  Ux_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdx);
  
  double
  Uy_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdy);
  
  double time, timed, timedd;
  double nu_htot = 0.;
  
  const distributed_vector& dg_coefficients;
  const distributed_vector& dg_coefficients_old;
  const distributed_vector& dg_coefficients_oldold;
  distributed_vector& incr_dg_coefficients;
  Q0& sol_onehalf;
  const Q1& Z;
  const Q0& slope_x;
  const Q0& slope_y;
  
private:

  std::array<double, 4> vel_rusanov_x, vel_rusanov_y, isdof_or_hanging, der_coeffs_x, der_coeffs_y;
  std::array<double, 2> grad_cell_h, grad_cell_Ux, grad_cell_Uy, grad_cell_ux, grad_cell_uy;
  
  const ordering& ordh;
  const ordering& ordUx;
  const ordering& ordUy;
  const ordering& ordh_nodal;
  const ordering& ordUx_nodal;
  const ordering& ordUy_nodal;
  const double& DELTAT;
  const double& epsilon;
  const bool& is_non_reflBC;
  const bool& is_bed_friction;
  const bool& is_stress_tensor;
  const double& grav;
  const double& density;
  const double& turbulence_coeff;
  const double& surface_pressure;
  const double& bed_friction_angle_rad;
  const double& fluid_viscosity;
  const double& yield_shear_stress;
  
};

double
signum (const double& x);


#endif
