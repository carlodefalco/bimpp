#ifndef TAYLOR_GALERKIN_H
#define TAYLOR_GALERKIN_H

#include <numeric> 
#include <bim_distributed_vector.h>
#include <tmesh.h>
#include <quad_operators.h>

 


class TG2_scheme  
{
  using Q1  = q1_vec<distributed_vector>;
  using Q0  = distributed_vector;
  
public:
  
  TG2_scheme(Q1& sol,
             Q1& sold,
             Q1& soldd,
             Q1& sold_rkc,
             Q1& soldd_rkc,
             Q1& sol_ini_rkc,
             Q1& incr,
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
             const Q1& Z,
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
             const std::vector<double>& extrema_vector,
             std::vector<std::array<double, 7> >& neig_state,
             const std::vector<double>& slope_x,
             const std::vector<double>& slope_y);
  
  TG2_scheme() = delete;
  
  ~TG2_scheme() = default;
  

  std::array<double,2>
  max_eigen (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  void
  compute_dt (tmesh::quadrant_iterator quadrant);
  
  void
  compute_dt_adaptive (tmesh::quadrant_iterator quadrant);

  void
  solve_non_lin(const int& kk);
  
  void
  first_step (tmesh::quadrant_iterator quadrant);

  void
  compute_nodal_anti_diffusive_fluxes (tmesh::quadrant_iterator quadrant);
  
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

  void
  communication_part(tmesh::quadrant_iterator quadrant, std::vector<MPI_Request>& reqs, int& count_req, int& shift);

  void
  communication_part(tmesh::quadrant_iterator quadrant, int& shift);
  
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
  std::array<double, 4> hwdof   = {0, 0, 0, 0};
  std::array<double, 4> hsdof   = {0, 0, 0, 0};
  std::array<double, 4> Uxwdof  = {0, 0, 0, 0};
  std::array<double, 4> Uywdof  = {0, 0, 0, 0};
  std::array<double, 4> Uxsdof  = {0, 0, 0, 0};
  std::array<double, 4> Uysdof  = {0, 0, 0, 0};
  std::array<double, 4> Z_node  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_hw_dof   = {0, 0, 0, 0};
  std::array<double, 4> P_minus_hw_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_hs_dof   = {0, 0, 0, 0};
  std::array<double, 4> P_minus_hs_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uxw_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uxw_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uyw_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uyw_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uxs_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uxs_dof = {0, 0, 0, 0};
  std::array<double, 4> P_plus_Uys_dof  = {0, 0, 0, 0};
  std::array<double, 4> P_minus_Uys_dof = {0, 0, 0, 0};
  
  std::array<double, 4> fluxx_hw_node    = {0, 0, 0, 0}, fluxy_hw_node    = {0, 0, 0, 0};
  std::array<double, 4> fluxx_hs_node    = {0, 0, 0, 0}, fluxy_hs_node    = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uxw_node   = {0, 0, 0, 0}, fluxy_Uxw_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uyw_node   = {0, 0, 0, 0}, fluxy_Uyw_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uxs_node   = {0, 0, 0, 0}, fluxy_Uxs_node   = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uys_node   = {0, 0, 0, 0}, fluxy_Uys_node   = {0, 0, 0, 0};
  
  
  
  // flux functions
  double
  hw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  hw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  hs_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  hs_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uxw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uxw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uyw_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uyw_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  Uxs_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uxs_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uys_flux_formula_x (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uys_flux_formula_y (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  



  // slope source terms
  double 
  src_slope_formula (const double& h, const double& S);

  double 
  src_slope_formula (const double& h, const double& S_x, const double& S_y, const int& kk);

  
  // source terms
  double
  hw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  hs_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uxs_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uys_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);

  double
  Uxw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);
  
  double
  Uyw_src_formula (const double& hw, const double& hs, const double& Uxw, const double& Uyw, const double& Uxs, const double& Uys);


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
  std::vector<std::array<double,4>>& incr_anti_diff;
  Q1& P_plus;
  Q1& P_minus;
  Q0& sol_onehalf;
  const Q1& Z;
  Q0& Z_onehalf;
  Q1& mass;
  
private:

  std::array<double, 4> vel_rusanov_x, vel_rusanov_y, isdof_or_hanging, der_coeffs_x, der_coeffs_y, der_coeffs_x_s, der_coeffs_y_s, D_U;
  std::array<double, 2> grad_cell_eta, grad_cell_hw, grad_cell_hs, grad_cell_Uxw, grad_cell_Uyw, grad_cell_Uxs, grad_cell_Uys, grad_cell_ux, grad_cell_uy, grad_cell_spec;
  
  const ordering& ordhw;
  const ordering& ordhs;
  const ordering& ordUxw;
  const ordering& ordUyw;
  const ordering& ordUxs;
  const ordering& ordUys;
  const double& DELTAT;
  const double& epsilon;
  const bool& is_non_reflBC;
  const double& grav;
  const double& erosion_coefficient;
  const double& density_w;
  const double& density_s;
  const double& turbulence_coeff;
  const bool& is_bed_friction;
  const double& bed_friction_angle_rad;
  const double& m_coeff;
  const double& terminal_velocity;
  const std::vector<double>& extrema_vector;
  std::vector<std::array<double, 7> >& neig_state;
  const std::vector<double>& slope_x;
  const std::vector<double>& slope_y;

  double w0, w1;

  // 5 arrays of storage as in Verwer's paper IMEX-RKCs,
  std::vector<double> b_vect, mu_tilde_vect, gamma_tilde_vect, v_vect, mu_vect;

  const double r_coeff = density_w/density_s;

  const double tol_incr = 1e-8;

  const double epsilon_IMEXRKC = 2./13.;

  const double tolerance_sign = 1.;

  const double regularization_parameter = 1e3; // has dimension of seconds, in this case the limit of the Bingham viscosity for small I_{2,D} exists finites
  
};




#endif
