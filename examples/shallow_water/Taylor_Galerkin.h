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
  
  TG2_scheme(const Q1& sol,
             const Q1& sold,
             const Q1& soldd,
             Q1& incr,
             Q0& sol_onehalf,
             const ordering& oh,
             const ordering& oUx, 
             const ordering& oUy, 
             const Q1& Z,
             const double& DELTAT,
             const double& h_min);
  
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
  flux_limiter(const double& Q_min, const double& Q_max, const double& Q_vertex, const double& Q_cell, const double& toll, double& phi_cell_Q);
  
  
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

  
  // source terms
  double
  h_src_formula (const double& h, const double& Ux, const double& Uy);
  
  double
  Ux_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdx);
  
  double
  Uy_src_formula (const double& h, const double& Ux, const double& Uy, const double& dZdy);
  
  double time, timed, timedd;
  double nu_htot = 0.;
  
  const Q1& sol;
  const Q1& sold;
  const Q1& soldd;
  Q1& incr;
  Q0& sol_onehalf;
  const Q1& Z;
  
private:
  static constexpr double grav = 9.81;
  
  std::array<double, 4> der_coeffs_x, der_coeffs_y;
  std::array<double, 4> vel_rusanov_x, vel_rusanov_y, isdof_or_hanging;
  double vel_rusanov_cell_x, vel_rusanov_cell_y;
  
  std::array<double, 2> grad_cell_h, grad_cell_Ux, grad_cell_Uy;
  
  
  const ordering& ordh;
  const ordering& ordUx;
  const ordering& ordUy;
  const double& DELTAT;
  const double& epsilon;
  
};




#endif
