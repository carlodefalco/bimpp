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
  
  TG2_scheme(Q1& sol_s,
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
             const double& erosion_coefficient);



  
  TG2_scheme() = delete;
  
  ~TG2_scheme() = default;
  
  
  void
  compute_dt (tmesh::quadrant_iterator quadrant);
  
  void
  compute_dt_adaptive (tmesh::quadrant_iterator quadrant);
  
  void
  first_step (tmesh::quadrant_iterator quadrant);

  void
  solve_non_lin (const int& kk);

  void
  compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant);

  void
  loop_step (const int& kk, const bool& isInitial);

  void
  loop_step_balance (tmesh::quadrant_iterator quadrant);

  void
  compute_stress_slope (tmesh::quadrant_iterator quadrant, const bool& isInitial);
  
  void
  second_step (tmesh::quadrant_iterator quadrant);
  
  void
  flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& mass_node, double& phi_cell_Q);

  void
  rkc(const int& j, const int& s, const int& kk);

  void
  set_dt (const double dt_);

  void
  set_old_dt (const double dt_);
  
  void
  set_times(const double& time, const double& time_old, const double& time_oldd);
  
  double
  get_dt ();
  
  
  double dt, dt_old;
  
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
  std::array<double, 4> hdof_s    = {0, 0, 0, 0};
  std::array<double, 4> Uxdof_s   = {0, 0, 0, 0};
  std::array<double, 4> Uydof_s   = {0, 0, 0, 0};
  std::array<double, 4> hdof_w    = {0, 0, 0, 0};
  std::array<double, 4> Uxdof_w   = {0, 0, 0, 0};
  std::array<double, 4> Uydof_w   = {0, 0, 0, 0};
  std::array<double, 4> Z_node  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_h_dof_s   = {0, 0, 0, 0};
  std::array<double, 4> P_minus_h_dof_s  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Ux_dof_s  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Ux_dof_s = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uy_dof_s  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uy_dof_s = {0, 0, 0, 0};
  std::array<double, 4> P_plus_h_dof_w   = {0, 0, 0, 0};
  std::array<double, 4> P_minus_h_dof_w  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Ux_dof_w  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Ux_dof_w = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uy_dof_w  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uy_dof_w = {0, 0, 0, 0};
  
  // std::array<double, 4> source_h_node  = {0, 0, 0, 0};
  // std::array<double, 4> source_Ux_node = {0, 0, 0, 0};
  // std::array<double, 4> source_Uy_node = {0, 0, 0, 0};
  
  std::array<double, 4> fluxx_h_node_s    = {0, 0, 0, 0}, fluxy_h_node_s    = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Ux_node_s   = {0, 0, 0, 0}, fluxy_Ux_node_s   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uy_node_s   = {0, 0, 0, 0}, fluxy_Uy_node_s   = {0, 0, 0, 0};

  std::array<double, 4> fluxx_h_node_w    = {0, 0, 0, 0}, fluxy_h_node_w    = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Ux_node_w   = {0, 0, 0, 0}, fluxy_Ux_node_w   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uy_node_w   = {0, 0, 0, 0}, fluxy_Uy_node_w   = {0, 0, 0, 0};
  
  
  std::array<double, 3> sigma_stress = {0., 0., 0.};
  
  
  // flux functions
  double
  h_flux_formula_x (const double& h, const double& Ux, const double& Uy);
  
  double
  h_flux_formula_y (const double& h, const double& Ux, const double& Uy);
  
  double
  Ux_flux_formula_x (const double& h, const double& Ux, const double& Uy, const double& h_);
  
  double
  Ux_flux_formula_y (const double& h, const double& Ux, const double& Uy);
  
  double
  Uy_flux_formula_x (const double& h, const double& Ux, const double& Uy);
  
  double
  Uy_flux_formula_y (const double& h, const double& Ux, const double& Uy, const double& h_);


  // stress functions
  double
  U_stress_formula (const double& h, const double& Ux, const double& Uy);
  


  std::array<double,3>
  compute_cell_stress (const double& Uxdof_0, const double& Uxdof_1, 
  const double& Uxdof_2, const double& Uxdof_3, 
  const double& Uydof_0, const double& Uydof_1, 
  const double& Uydof_2, const double& Uydof_3);

  std::array<double,6>
  compute_cell_def_grad (const double& Uxdof_0, const double& Uxdof_1, 
  const double& Uxdof_2, const double& Uxdof_3, 
  const double& Uydof_0, const double& Uydof_1, 
  const double& Uydof_2, const double& Uydof_3);

  // slope source terms
  double 
  src_slope_formula (const double& h, const double& S);

  double 
  src_slope_formula (const double& h, const double& S_x, const double& S_y, const int& kk);

  
  // source terms
  double
  h_src_formula (const double& h, const double& Ux, const double& Uy, const double& n_poro);
  
  double
  friction_Ux (const double& h, const double& Ux, const double& Uy, const double& n_poro);
  
  double
  friction_Uy (const double& h, const double& Ux, const double& Uy, const double& n_poro);

  double
  Ux_src_formula_w (const double& h, const double& U_wx, const double& U_sx, const double& U_wy, const double& U_sy, const double& n_poro);

  double
  Ux_src_formula_s (const double& h, const double& U_wx, const double& U_sx, const double& U_wy, const double& U_sy, const double& n_poro);

  double
  Uy_src_formula_w (const double& h, const double& U_wx, const double& U_sx, const double& U_wy, const double& U_sy, const double& n_poro);

  double
  Uy_src_formula_s (const double& h, const double& U_wx, const double& U_sx, const double& U_wy, const double& U_sy, const double& n_poro);

  double
  int_term_x (const double& h, const double& U_wy, const double& U_sy, const double& n_poro, const double& dens_a);

  double
  int_term_y (const double& h, const double& U_wy, const double& U_sy, const double& n_poro, const double& dens_a);

  void
  prepare_IMEXRKC_coefficients (const int& s);
  
  double time, timed, timedd;
  double nu_htot = 0.;
  
  Q1& sol_s;
  Q1& sol_w;
  Q1& sold_s;
  Q1& sold_w;
  Q1& soldd_s;
  Q1& soldd_w;
  Q1& sold_rkc_s;
  Q1& sold_rkc_w;
  Q1& soldd_rkc_s;
  Q1& soldd_rkc_w;
  Q1& sol_ini_rkc_s;
  Q1& sol_ini_rkc_w;
  Q1& incr_s;
  Q1& incr_w;
  Q1& incr_initial_source_s;
  Q1& incr_initial_source_w;
  Q1& incr_source_s;
  Q1& incr_source_w;
  std::vector<std::array<double,4>>& incr_anti_diff_s;
  std::vector<std::array<double,4>>& incr_anti_diff_w;
  Q1& P_plus_s;
  Q1& P_plus_w;
  Q1& P_minus_s;
  Q1& P_minus_w;
  Q1& spec_radius_nodal_s;
  Q1& spec_radius_nodal_w;
  Q0& sol_onehalf_s;
  Q0& sol_onehalf_w;
  Q1& porosity;
  const Q1& Z;
  Q1& Newton_it;
  Q1& slope_x_node;
  Q1& slope_y_node;
  Q1& stress_initial_step_s;
  Q1& stress_initial_step_w;  
  Q1& stress_step_s;
  Q1& stress_step_w;
  const Q0& slope_x;
  const Q0& slope_y;
  Q1& mass;
  
private:

  std::array<double, 4> vel_rusanov_x, vel_rusanov_y, isdof_or_hanging, der_coeffs_x, der_coeffs_y, der_coeffs_x_s, der_coeffs_y_s, D_U;
  std::array<double, 2> grad_cell_eta_s, grad_cell_h_s, grad_cell_Ux_s, grad_cell_Uy_s, grad_cell_eta_w, grad_cell_h_w, grad_cell_Ux_w, grad_cell_Uy_w, grad_cell_ux, grad_cell_uy, grad_cell_spec;
  
  const ordering& ordh;
  const ordering& ordUx;
  const ordering& ordUy;
  const double& DELTAT;
  const double& epsilon;
  const bool& is_non_reflBC;
  const bool& is_bed_friction;
  const bool& is_stress_tensor; 
  const bool& is_erosion;
  const double& grav;
  const double& density;
  const double& density_s;
  const double& density_w;
  const double& turbulence_coeff;
  const double& surface_pressure;
  const double& bed_friction_angle_rad;
  const double& fluid_viscosity;
  const double& yield_shear_stress;
  const double& erosion_coefficient;

  double w0, w1;

  // 5 arrays of storage as in Verwer's paper IMEX-RKCs,
  std::vector<double> b_vect, mu_tilde_vect, gamma_tilde_vect, v_vect, mu_vect;

  const double tol_incr = 1e-8;

  const double epsilon_IMEXRKC = 2./13.;

  const double tolerance_sign = 1e-2;

  const double regularization_parameter = 1e3; // has dimension of seconds, in this case the limit of the Bingham viscosity for small I_{2,D} exists finites
  
};




#endif
