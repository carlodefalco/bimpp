#ifndef TAYLOR_GALERKIN_H
#  define TAYLOR_GALERKIN_H

#include <numeric>
#include <bim_distributed_vector.h>
#include <tmesh.h>
#include <quad_operators.h>

 


class TG2_scheme
{
  using Q1  = q1_vec<distributed_vector>;
  using Q0  = std::vector<double>;
  
public:
  
  TG2_scheme(const Q1& state,
             Q0& state_onehalf,
             const ordering& oh,
             const ordering& oUx,
             const ordering& oUy,
             const Q1& elevation,
             const double& DELTAT);
  
  TG2_scheme() = delete;
  
  ~TG2_scheme() = default;
  
  void
  set_quadrant (tmesh::quadrant_iterator quadrant);
  
  void
  step_function (Q1& increment, tmesh::quadrant_iterator quadrant);

  
  
  void
  compute_dt (tmesh::quadrant_iterator quadrant);
  
  void
  first_step (tmesh::quadrant_iterator quadrant, const double& deltat);
  
  void
  second_step (tmesh::quadrant_iterator quadrant, Q1& increment);
  
  
  void
  set_dt (const double dt_);
  
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
  std::array<double, 4> hdof  = {0, 0, 0, 0};
  std::array<double, 4> Uxdof = {0, 0, 0, 0};
  std::array<double, 4> Uydof = {0, 0, 0, 0};
  
  std::array<double, 4> source_h_node  = {0, 0, 0, 0};
  std::array<double, 4> source_Ux_node = {0, 0, 0, 0};
  std::array<double, 4> source_Uy_node = {0, 0, 0, 0};
  
  std::array<double, 4> fluxx_h_node  = {0, 0, 0, 0}, fluxy_h_node  = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Ux_node = {0, 0, 0, 0}, fluxy_Ux_node = {0, 0, 0, 0};
  std::array<double, 4> fluxx_Uy_node = {0, 0, 0, 0}, fluxy_Uy_node = {0, 0, 0, 0};
  
  
  
  
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
  Ux_src_formula (const double& h, const double& Ux, const double& Uy);
  
  double
  Uy_src_formula (const double& h, const double& Ux, const double& Uy);

  
  Q0& sol_onehalf;
  
private:
  static constexpr double grav = 9.81;
  static constexpr double epsilon = 1.e-6;
  
  std::array<double, 4> der_coeffs_x, der_coeffs_y;
  std::array<double, 4> vel_rusanov_x, vel_rusanov_y, isdof_or_hanging;
  
  int gn_elements;
  
  const Q1& state_vector;
  
  const ordering& ordh;
  const ordering& ordUx;
  const ordering& ordUy;
  const Q1& Z;
  const double& DELTAT;
  
};




#endif
