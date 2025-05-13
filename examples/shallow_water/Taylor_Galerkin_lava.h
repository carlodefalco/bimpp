#ifndef TAYLOR_GALERKIN_H
#define TAYLOR_GALERKIN_H

#include <numeric> 
#include <bim_distributed_vector.h>
#include <tmesh.h>
#include <quad_operators.h>


#define SET_COEFFICIENTS 3

class TG2_scheme  
{
  using Q1  = q1_vec<distributed_vector>;
  using Q0  = distributed_vector;
  
public:
  
  TG2_scheme(Q1& sol,
             Q1& sold,
             Q1& soldd,
             Q1& sol_2,
             Q1& incr,
             Q1& incr_second,
             std::vector<std::array<double,4>>& incr_anti_diff,
             Q1& P_plus,
             Q1& P_minus,
             Q0& sol_onehalf,
             Q1& mass,
             const ordering& oh,
             const ordering& oUx, 
             const ordering& oUy, 
             const ordering& oTh, 
             const Q1& Z,
             Q0& Z_onehalf,
             const double& DELTAT,
             const double& h_min,
             const bool& is_non_reflBC,
             const bool& is_isothermal,
             const bool& is_limiter,
             const double& grav,
             const double& nu_ref,
             const double& T_ref,
             const double& b_exp_coeff,
             const double& saturation_coeff,
             const double& density,
	           const double& T_env,
             const double& specific_heat_pressure,
             const double& convective_coeff,
             const double& x_v,
             const double& y_v,
             const double& Q_vent, 
             const double& T_vent,
             const double& sigma_vent);
  
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
  loop_step (const int& kk, const bool& isInitial);

  void
  compute_stress_slope (tmesh::quadrant_iterator quadrant, const bool& isInitial);
  
  void
  second_step (tmesh::quadrant_iterator quadrant);
  
  void
  flux_limiter(const double& Q_min, const double& Q_max, const double& Q_dof, const double& P_plus_Q, const double& P_minus_Q, const double& flux_on_the_node, const double& mass_node, double& phi_cell_Q);

  void
  rkc(const int& j, const int& s, const int& kk);

  double
  Ux_jac_source(const double& h, const double& Ux, const double& Uy);

  double
  Uy_jac_source(const double& h, const double& Ux, const double& Uy);

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
  
  
  double dt, dt_old, dt_22, dt_33, dt_21, dt_31, dt_32, dt_expl_21, dt_expl_32, b_1, b_2, b_3, b_expl_1, b_expl_2, b_expl_3;
  
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
  std::array<double, 4> Uxdof   = {0, 0, 0, 0};
  std::array<double, 4> Uydof   = {0, 0, 0, 0};
  std::array<double, 4> Thdof   = {0, 0, 0, 0};
  std::array<double, 4> Z_node  = {0, 0, 0, 0};
  std::array<double, 4> Z_node_nei  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_h_dof   = {0, 0, 0, 0};
  std::array<double, 4> P_minus_h_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Ux_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Ux_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uy_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uy_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Th_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Th_dof = {0, 0, 0, 0};
  
  // std::array<double, 4> source_h_node  = {0, 0, 0, 0};
  // std::array<double, 4> source_Ux_node = {0, 0, 0, 0};
  // std::array<double, 4> source_Uy_node = {0, 0, 0, 0};
  
  std::array<double, 4> fluxx_h_node    = {0, 0, 0, 0}, fluxy_h_node    = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Ux_node   = {0, 0, 0, 0}, fluxy_Ux_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uy_node   = {0, 0, 0, 0}, fluxy_Uy_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Th_node   = {0, 0, 0, 0}, fluxy_Th_node   = {0, 0, 0, 0};
  
  
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

  double
  Th_flux_formula_x (const double& h, const double& Ux, const double& Uy, const double& Th);
  
  double
  Th_flux_formula_y (const double& h, const double& Ux, const double& Uy, const double& Th);


  // slope source terms
  double 
  src_slope_formula (const double& h, const double& S);

  void
  solve_non_lin(const int& kk);

  void
  compute_updated_sol(tmesh::quadrant_iterator quadrant);

  void
  compute_updated_sol(const int& kk);

  
  // source terms
  double
  h_src_formula (const double& h, const double& Ux, const double& Uy);
  
  double
  Ux_src_formula (const double& h, const double& Ux, const double& Uy, const double& Th);
  
  double
  Uy_src_formula (const double& h, const double& Ux, const double& Uy, const double& Th);

  double
  Th_src_formula (const double& h, const double& Ux, const double& Uy, const double& Th, const double& T_enva);

  void
  prepare_IMEXRKC_coefficients (const int& s);

  double
  signum (const double& x);
  
  double time, timed, timedd;
  double nu_htot = 0.;
  double Fr = 0.;
  double g_coeff = 0.;
  
  Q1& sol;
  Q1& sold;
  Q1& soldd;
  Q1& sol_2;
  Q1& incr;
  Q1& incr_second;
  std::vector<std::array<double,4>>& incr_anti_diff;
  Q1& P_plus;
  Q1& P_minus;
  Q0& sol_onehalf;
  const Q1& Z;
  Q0& Z_onehalf;
  Q1& mass;
  
private:

  std::array<double, 4> vel_rusanov_x, vel_rusanov_y, isdof_or_hanging, der_coeffs_x, der_coeffs_y, der_coeffs_x_s, der_coeffs_y_s, D_U;
  std::array<double, 2> grad_cell_Z, grad_cell_eta, grad_cell_h, grad_cell_Ux, grad_cell_Uy, grad_cell_Th, grad_cell_ux, grad_cell_uy, grad_cell_spec;
  
  const ordering& ordh;
  const ordering& ordUx;
  const ordering& ordUy;
  const ordering& ordTh;
  const double& DELTAT;
  const double& epsilon;
  const bool& is_non_reflBC;
  const bool& is_isothermal;
  const bool& is_limiter;
  const double& grav;
  const double& nu_ref;
  const double& T_ref;
  const double& b_exp_coeff;
  const double& saturation_coeff;
  const double& density;
  const double& T_env;                     
  const double& specific_heat_pressure;
  const double& convective_coeff;
  const double& x_v;
  const double& y_v;
  const double& Q_vent;
  const double& T_vent;
  const double& sigma_vent;

  double w0, w1;

  // 5 arrays of storage as in Verwer's paper IMEX-RKCs,
  std::vector<double> b_vect, mu_tilde_vect, gamma_tilde_vect, v_vect, mu_vect;

  const double tol_incr = 1e-8;

  const double epsilon_IMEXRKC = 2./13.;

  const double tolerance_sign = .1;

  const double regularization_parameter = 1e3; // has dimension of seconds, in this case the limit of the Bingham viscosity for small I_{2,D} exists finites
  
};




#endif
