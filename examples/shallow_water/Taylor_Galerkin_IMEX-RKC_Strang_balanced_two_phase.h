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
  
  TG2_scheme(Q1& sol,
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
             const double& terminal_velocity);
  
  TG2_scheme() = delete;
  
  ~TG2_scheme() = default;
  
  
  void
  compute_dt (tmesh::quadrant_iterator quadrant);
  
  void
  compute_dt_adaptive (tmesh::quadrant_iterator quadrant);
  
  void
  first_step (tmesh::quadrant_iterator quadrant);

  void
  compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant);

  void
  solve_non_lin(const int& kk);

  void
  low_order_sol(const int& kk);

  void
  loop_step (const int& kk, const bool& isInitial);

  void
  compute_stress_slope (tmesh::quadrant_iterator quadrant, const bool& isInitial);
  
  void
  second_step (tmesh::quadrant_iterator quadrant);
  
  void
  flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& mass_node, double& phi_cell_Q);

  void
  rkc(const int& j, const int& s, const int& kk);



  /*
  void
  mu_fun(const int& s);

  void
  v_fun(const int& s);

  void
  gamma_tilde_fun(const int& s);

  void
  T_fun_second(const int& s);
  
  void
  T_fun_prime(const int& s);

  void
  T_fun(const int& s);

  void
  w_fun_0(const int& s);

  void
  w_fun_1(const int& s);

  void
  b_fun(const int& s);

  void
  mu_fun_tilde (const int& s);
  */

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
  std::array<double, 4> etadof  = {0, 0, 0, 0};
  std::array<double, 4> hdof    = {0, 0, 0, 0};
  std::array<double, 4> ndof    = {0, 0, 0, 0};
  std::array<double, 4> Uxwdof  = {0, 0, 0, 0};
  std::array<double, 4> Uywdof  = {0, 0, 0, 0};
  std::array<double, 4> Uxsdof  = {0, 0, 0, 0};
  std::array<double, 4> Uysdof  = {0, 0, 0, 0};
  std::array<double, 4> Z_node  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_h_dof    = {0, 0, 0, 0};
  std::array<double, 4> P_minus_h_dof   = {0, 0, 0, 0};
  std::array<double, 4> P_plus_n_dof    = {0, 0, 0, 0};
  std::array<double, 4> P_minus_n_dof   = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uxw_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uxw_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uyw_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uyw_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uxs_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uxs_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uys_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uys_dof = {0, 0, 0, 0};
  
  // std::array<double, 4> source_h_node  = {0, 0, 0, 0};
  // std::array<double, 4> source_Ux_node = {0, 0, 0, 0};
  // std::array<double, 4> source_Uy_node = {0, 0, 0, 0};
  
  std::array<double, 4> fluxx_h_node     = {0, 0, 0, 0}, fluxy_h_node     = {0, 0, 0, 0};
  std::array<double, 4> fluxx_n_node     = {0, 0, 0, 0}, fluxy_n_node     = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uwx_node   = {0, 0, 0, 0}, fluxy_Uwx_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uwy_node   = {0, 0, 0, 0}, fluxy_Uwy_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Usx_node   = {0, 0, 0, 0}, fluxy_Usx_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Usy_node   = {0, 0, 0, 0}, fluxy_Usy_node   = {0, 0, 0, 0};
  
  
  std::array<double, 3> sigma_stress = {0., 0., 0.};
  
  
  // flux functions
  double
  h_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  h_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  n_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  n_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uwx_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uwx_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uwy_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uwy_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  Usx_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Usx_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Usy_flux_formula_x (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Usy_flux_formula_y (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);


  // stress functions
  double
  U_stress_formula (const double& h, const double& Ux, const double& Uy);


  // max eigenvalues func
  std::array<double,2>
  max_eigen (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  


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
  h_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uxs_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uys_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  Uxw_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uyw_src_formula (const double& h, const double& n, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);


  void
  prepare_IMEXRKC_coefficients (const int& s);

  double
  signum (const double& x);
  
  double time, timed, timedd;
  double nu_htot = 0.;
  
  Q1& sol;
  Q1& sold;
  Q1& soldd;
  Q1& sold_rkc;
  Q1& soldd_rkc;
  Q1& sol_ini_rkc;
  Q1& incr;
  Q1& incr_initial_source;
  Q1& incr_source;
  std::vector<std::array<double,4>>& incr_anti_diff;
  Q1& P_plus;
  Q1& P_minus;
  Q1& spec_radius_nodal;
  Q0& sol_onehalf;
  const Q1& Z;
  Q0& Z_onehalf;
  Q1& Newton_it;
  Q1& stress_initial_step; 
  Q1& stress_step;
  Q1& mass;
  
private:

  std::array<double, 4> vel_rusanov_x, vel_rusanov_y, isdof_or_hanging, der_coeffs_x, der_coeffs_y, der_coeffs_x_s, der_coeffs_y_s, D_U;
  std::array<double, 2> grad_cell_eta, grad_cell_h, grad_cell_n, grad_cell_Uxw, grad_cell_Uyw, grad_cell_Uxs, grad_cell_Uys, grad_cell_ux, grad_cell_uy, grad_cell_spec;
  
  const ordering& ordh;
  const ordering& ordn;
  const ordering& ordUxw;
  const ordering& ordUyw;
  const ordering& ordUxs;
  const ordering& ordUys;
  const double& DELTAT;
  const double& epsilon;
  const bool& is_non_reflBC;
  const bool& is_stress_tensor; 
  const double& grav;
  const double& erosion_coefficient;
  const double& density_w;
  const double& density_s;
  const double& turbulence_coeff;
  const double& bed_friction_angle_rad;
  const double& fluid_viscosity;
  const double& yield_shear_stress;
  const double& m_coeff;
  const double& terminal_velocity;

  double w0, w1;

  // 5 arrays of storage as in Verwer's paper IMEX-RKCs,
  std::vector<double> b_vect, mu_tilde_vect, gamma_tilde_vect, v_vect, mu_vect;

  const double tol_incr = 1e-8;

  const double epsilon_IMEXRKC = 2./13.;

  const double tolerance_sign = 1e-2;

  const double regularization_parameter = 1e3; // has dimension of seconds, in this case the limit of the Bingham viscosity for small I_{2,D} exists finites
  
};




#endif
