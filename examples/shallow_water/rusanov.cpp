/*
  Copyright (C) 2020 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <octave_file_io.h>

#include <bim_distributed_vector.h>
#include <bim_timing.h>
#include <tmesh.h>
#include <quad_operators.h>

static constexpr char LOADFILENAME_1[255] = "orography.octbin.gz";
static constexpr char LOADFILENAME_2[255] = "basin.octbin.gz";
static constexpr char SAVEFILENAME_1[255] = "orography_tmsh";
static constexpr char VARNAME_1[255] = "dem";
static constexpr char VARNAME_2[255] = "basin";

// properties of the input dem
static constexpr double res = 5; // it is also the minimum resolution of the bim element
static constexpr double Nx = 5;//1998;
static constexpr double Ny = 5;//1829;

// scale bim domain st have only multiple of res has bim elements
// L, H are fictitious
// assert that std::ceil(std::log2(Nx-1)=std::ceil(std::log2(Ny-1)
//static constexpr double max_num_uniform_refinement = std::ceil(std::log2(Nx-1));

static constexpr double L = res*(Nx-1);
static constexpr double H = res*(Ny-1);
//static constexpr double scale_factor_x = std::pow(2., max_num_uniform_refinement)/(Nx-1);
//static constexpr double scale_factor_y = std::pow(2., max_num_uniform_refinement)/(Ny-1);
static std::vector<double>   dem;
static std::vector<double>   basin_mask;
static constexpr int NUM_REFINEMENTS = 6; // 3, to get minimum refinement, i.e. bim element reolution equal to res on the whole domain, put std::pow(2., std::ceil(std::log2(Ny-1))=std::pow(2., std::ceil(std::log2(Nx-1))
static constexpr int NUM_TREFINEMENTS = 1; // 10


static constexpr double SAVEDT  = 1e-3;
static constexpr double DELTAT =  1e-3;
static constexpr double REDCDT =  1;
static constexpr double T      =  5;





// Connectivity of local element
constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {0,  0,
   0,  H,
   L,  0,
   L,  H};

const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {  1,    3,    4,    2,    1 };

// Refinement rule
static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
hanging_refinement (tmesh::quadrant_iterator quadrant)
{
  double x_minus, x_plus, y_minus, y_plus;
  x_minus = quadrant->p (0, 0);
  x_plus  = quadrant->p (0, 1);
  y_minus = quadrant->p (1, 0);
  y_plus  = quadrant->p (1, 2);
  
  const auto marker = (x_minus+x_plus)/2<L/2;//(x_minus+x_plus)/2<L/2 && (y_minus+y_plus)/2>H/2 ? 1 : 0;
  return marker; }


static int
raster_2_vector(const double& i_x,
                const double& i_y)
{
  static int ii;
  ii = i_y + Ny*i_x;
  
  return(ii);
}

static std::array<int,3>
global_coord_2_raster(const double& x,
                      const double& y)
{
  static double i_x;
  static double i_y;
  static int ii;
  
  
  // nearest neighbor
  i_x = std::round(x / res);
  i_y = - std::round(y / res) + (Ny-1);
  
  ii = raster_2_vector(i_x,i_y);
  
  return(std::array<int,3>{{ ii,int(i_x),int(i_y) }});
}

static double
refine_function (tmesh::quadrant_iterator quadrant)
{
  
  static double xm;
  static double ym;
  
  static std::array<int,3> ids;
  
  static double x_minus;
  static double x_plus;
  static double y_minus;
  static double y_plus;
  
  static double grad_x;
  static double grad_y;
  static double grad;
  static double grad_dem;
  
  static double rel_error;
  static double N_el, N_x, N_y;
  
  static double res_x, res_y;
  static double Dx_bim, Dy_bim;
  
  static double const toll = 0.2;//0.002; //0.2
  
  x_minus = quadrant->p (0, 0);
  x_plus  = quadrant->p (0, 1);
  y_minus = quadrant->p (1, 0);
  y_plus  = quadrant->p (1, 2);
  
  const auto ids_minus_minus = global_coord_2_raster(x_minus,y_minus);
  const auto ids_minus_plus  = global_coord_2_raster(x_minus,y_plus);
  const auto ids_plus_minus  = global_coord_2_raster(x_plus,y_minus);
  const auto ids_plus_plus   = global_coord_2_raster(x_plus,y_plus);
  
  const auto & i_x_m = ids_minus_minus[1];
  const auto & i_x_p = ids_plus_minus[1];
  const auto & i_y_p = ids_minus_plus[2];
  const auto & i_y_m = ids_minus_minus[2];
  
  Dx_bim = (x_plus - x_minus);
  Dy_bim = (y_plus - y_minus);
  
  // gradient on the bimpp element
  grad_x = std::abs(dem[raster_2_vector(i_x_m,i_y_p)] - dem[raster_2_vector(i_x_p,i_y_p)]);
  grad_x /= Dx_bim;
  
  grad_y = std::abs(dem[raster_2_vector(i_x_m,i_y_p)] - dem[raster_2_vector(i_x_m,i_y_m)]);
  grad_y /= Dy_bim;
  
  grad = std::sqrt(std::pow(grad_x,2.) + std::pow(grad_y,2.));
  
  N_x = i_x_p-i_x_m;
  N_y = i_y_m-i_y_p;
  N_el = N_x*N_y;
  
  res_x = Dx_bim/N_x;
  res_y = Dy_bim/N_y;
  
  grad_dem = 0;
  rel_error = 0;
  for (int i_x = i_x_m; i_x < i_x_p; i_x++)
  {
    for (int i_y = i_y_p; i_y < i_y_m; i_y++)
    {
      const auto k  = raster_2_vector(i_x,i_y);
      const auto kx = raster_2_vector(i_x+1,i_y);
      const auto ky = raster_2_vector(i_x,i_y+1);
      
      const auto & h_center = dem[k];
      const auto & h_plus_x = dem[kx];
      const auto & h_plus_y = dem[ky];
      
      grad_x = std::abs(h_center - h_plus_x)/res_x;
      grad_y = std::abs(h_center - h_plus_y)/res_y;
      
      grad_dem = std::sqrt(std::pow(grad_x,2.) + std::pow(grad_y,2.));
      
      rel_error += std::pow(grad-grad_dem,2.);
      
    }
  }
  
//  rel_error = std::sqrt(N_el * rel_error * res);
  rel_error = std::sqrt(rel_error/(N_el==0 ? 1. : N_el));
  
  const auto basin_check = basin_mask[ids_minus_minus[0]] +
  basin_mask[ids_minus_plus[0]] + basin_mask[ids_plus_minus[0]] +
  basin_mask[ids_plus_plus[0]];
  
  
  return((basin_check>0 && ((rel_error>toll)*((x_plus-x_minus)>res))) ? 1 : 0);
  
}



static int
coarsen_function (tmesh::quadrant_iterator quadrant)
{
  static double xm;
  static double ym;
  
  static std::array<int,3> ids;
  
  static double x_minus;
  static double x_plus;
  static double y_minus;
  static double y_plus;
  
  static double grad_x;
  static double grad_y;
  static double grad;
  static double grad_dem;
  
  static double rel_error;
  
  x_minus = quadrant->p (0, 0);
  x_plus  = quadrant->p (0, 1);
  y_minus = quadrant->p (1, 0);
  y_plus  = quadrant->p (1, 2);
  
  const auto & ids_minus_minus = global_coord_2_raster(x_minus,y_minus);
  const auto & ids_minus_plus  = global_coord_2_raster(x_minus,y_plus);
  const auto & ids_plus_minus  = global_coord_2_raster(x_plus,y_minus);
  const auto & ids_plus_plus   = global_coord_2_raster(x_plus,y_plus);
  
  
  return(!basin_mask[ids_minus_minus[0]] &&
         !basin_mask[ids_minus_plus[0]] &&
         !basin_mask[ids_plus_minus[0]] &&
         !basin_mask[ids_plus_plus[0]] ? 1 : 0);
}




using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = std::vector<double>;         // Typedef for local q_0 vector

//double h0_fun (const double& xx, const double& yy)  { return std::max (0., (8. - std::sin (M_PI * xx / 2. / 400.) - dem[global_coord_2_raster(xx,yy)[0]])); }
double h0_fun (const double& xx, const double& yy)  {
  if (xx>L*3./10. && xx<7./10.*L)
  {
//    std::cout << "aa" << std::endl;
    return 2;
  }
  
//  if (xx>L/4 && xx<3/4*L && yy>H/4 && yy <3/4*H)
//  {
//    return 2;
//  }
  return 0;
}
double Ux0_fun (double xx, double yy) { return 0.; }
double Uy0_fun (double xx, double yy) { return 0.; }

// Assemble vector from mesh.
// FIXME  the following two functions are copied over from
// "quad_operators.cpp" as they were not exported in an header,
// should find better way to avoid code duplication
void
assemble_vector (tmesh::quadrant_iterator& quadrant,
                 const std::array<double, 4>& locrhs,
                 Q1& rhs,
                 const ordering& ord = default_ord)
{
  
  std::vector<unsigned int> rows;
  rows.reserve (2);
  int i, r;
  
  for (i = 0; i < 4; ++i)
  {
    rows.clear ();
    
    if (! quadrant->is_hanging (i))
      rows.push_back (quadrant->gt (i));
    else
    {
      rows.push_back (quadrant->gparent (0, i));
      rows.push_back (quadrant->gparent (1, i));
    }
    
    for (r = 0; r < rows.size (); ++r)
    rhs[ord (rows[r])] +=
    locrhs[i] / rows.size ();
  }
}

using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = std::vector<double>;         // Typedef for local q_0 vector


class
TG2_rusanov
{

private:
  static constexpr double grav = 9.81;
  
  std::array<double, 4> dZdx  = {0, 0, 0, 0};
  std::array<double, 4> dZdy  = {0, 0, 0, 0};
  
  std::array<double, 4> der_coeffs_x, der_coeffs_y,
  vel_rusanov_x, vel_rusanov_y, isdof_or_hanging;
  
  const Q1& state_vector;
  const ordering& ordh;
  const ordering& ordUx;
  const ordering& ordUy;
  const Q1& Z;

public:
//  using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
//  using Q0  = std::vector<double>;         // Typedef for local q_0 vector

  ~TG2_rusanov () = default;
  TG2_rusanov(const Q1& state,
              const ordering& oh,
              const ordering& oUx,
              const ordering& oUy,
              const Q1& elevation)
  : state_vector(state), ordh(oh), ordUx(oUx), ordUy(oUy), Z(elevation)
  {  };
  
  void
  set_quadrant (tmesh::quadrant_iterator quadrant)
  {
    for (int ii = 0; ii < 4; ++ii) {
      xn[ii] = quadrant->p(0, ii);
      yn[ii] = quadrant->p(1, ii);
    }
    
    Dx = xn[1]-xn[0];
    Dy = yn[2]-yn[0];
    area = Dx * Dy;
    
    
    for (int ii = 0; ii < 4; ++ii){
      if (! quadrant->is_hanging (ii)){
        hdof[ii] = state_vector [ordh (quadrant->gt (ii))];
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
    
    
    
    
    for (int ii = 0; ii < 4; ++ii){
      const auto& hpoint = hdof[ii];
      const auto celerity = std::sqrt(grav*hpoint);
      vel_rusanov_x[ii] = hpoint>0 ? std::abs(Uxdof[ii]/hpoint)+celerity : 0.;
      vel_rusanov_y[ii] = hpoint>0 ? std::abs(Uydof[ii]/hpoint)+celerity : 0.;
      
      const auto dtoptx = hpoint>0 ? Dx/vel_rusanov_x[ii] : DELTAT;
      const auto dtopty = hpoint>0 ? Dy/vel_rusanov_y[ii] : DELTAT;
      const auto dtopt = dtoptx > dtopty ? dtopty : dtoptx;
      if (dt > dtopt) set_dt (dtopt);
    }
    
    
    /*
    const double DxDyby4 = .25 * area;
    wq = {DxDyby4, DxDyby4, DxDyby4, DxDyby4};
    
    shgx = {
      {-1./Dx, -1./Dx,  0,      0},
      {1./Dx,   1./Dx,  0,      0},
      {0,       0,     -1./Dx, -1./Dx},
      {0,       0,      1./Dx,  1./Dx}
    };

    shgy = {
      {-1./Dy,  0,     -1./Dy,  0},
      { 0,     -1./Dy,  0,     -1./Dy},
      { 1./Dy,  0,      1./Dy,  0},
      { 0,      1./Dy,  0,      1./Dy}
    }
    
    dZdx.fill (0.0);
    dZdy.fill (0.0);
    for (int kk = 0; kk < 4; ++kk) {
      for (int ii = 0; ii < 4; ++ii){
        if (! quadrant->is_hanging (ii)){
          dZdx[kk] += shgx[ii][kk] * Z[quadrant->gt (ii)];
          dZdy[kk] += shgy[ii][kk] * Z[quadrant->gt (ii)];
        } else {
          dZdx[ii] = shgx[ii][kk] * .5 * (Z[quadrant->gparent(0,ii)] +
                                          Z[quadrant->gparent(1,ii)]);
          dZdy[ii] = shgy[ii][kk] * .5 * (Z[quadrant->gparent(0,ii)] +
                                          Z[quadrant->gparent(1,ii)]);
        }
      }
    }*/

  }
  
  void
  step_function (Q1& increment, tmesh::quadrant_iterator quadrant) {
    // in questa fz bisogna riempire incr
    int ii;
    
    double F_star_h_x = 0, F_star_h_y = 0,
    F_star_Ux_x = 0, F_star_Ux_y = 0,
    F_star_Uy_x = 0, F_star_Uy_y = 0;
    
    for (ii = 0; ii < 4; ++ii)
    {
      F_star_h_x += h_flux_formula_x  (hdof[ii], Uxdof[ii], Uydof[ii]);
      F_star_h_y += h_flux_formula_y  (hdof[ii], Uxdof[ii], Uydof[ii]);
      
      F_star_Ux_x += Ux_flux_formula_x (hdof[ii], Uxdof[ii], Uydof[ii]);
      F_star_Ux_y += Ux_flux_formula_y (hdof[ii], Uxdof[ii], Uydof[ii]);
      
      F_star_Uy_x += Uy_flux_formula_x (hdof[ii], Uxdof[ii], Uydof[ii]);
      F_star_Uy_y += Uy_flux_formula_y (hdof[ii], Uxdof[ii], Uydof[ii]);
    }
    F_star_h_x  /= 4.;
    F_star_h_y  /= 4.;
    F_star_Ux_x /= 4.;
    F_star_Ux_y /= 4.;
    F_star_Uy_x /= 4.;
    F_star_Uy_y /= 4.;
    
    
    // Lax-Friedrichs correction term, with Rusanov treatment
    F_star_h_x -= .5 * ( (vel_rusanov_y[3]*hdof[3] - vel_rusanov_y[2]*hdof[2]) +
                         (vel_rusanov_y[1]*hdof[1] - vel_rusanov_y[0]*hdof[0]) );
    F_star_h_y -= .5 * ( (vel_rusanov_x[2]*hdof[2] - vel_rusanov_x[0]*hdof[0]) +
                         (vel_rusanov_x[3]*hdof[3] - vel_rusanov_x[1]*hdof[1]) );
    
    F_star_Ux_x -= .5 * ( (vel_rusanov_y[3]*Uxdof[3] - vel_rusanov_y[2]*Uxdof[2]) +
                          (vel_rusanov_y[1]*Uxdof[1] - vel_rusanov_y[0]*Uxdof[0]) );
    F_star_Ux_y -= .5 * ( (vel_rusanov_x[2]*Uxdof[2] - vel_rusanov_x[0]*Uxdof[0]) +
                          (vel_rusanov_x[3]*Uxdof[3] - vel_rusanov_x[1]*Uxdof[1]) );
    
    F_star_Uy_x -= .5 * ( (vel_rusanov_y[3]*Uydof[3] - vel_rusanov_y[2]*Uydof[2]) +
                          (vel_rusanov_y[1]*Uydof[1] - vel_rusanov_y[0]*Uydof[0]) );
    F_star_Uy_y -= .5 * ( (vel_rusanov_x[2]*Uydof[2] - vel_rusanov_x[0]*Uydof[0]) +
                          (vel_rusanov_x[3]*Uydof[3] - vel_rusanov_x[1]*Uydof[1]) );
    
    
    
    for (int ii = 0; ii < 4; ++ii){
      
      const auto h_  = der_coeffs_x[ii] * F_star_h_x  + der_coeffs_y[ii] * F_star_h_y;
      const auto Ux_ = der_coeffs_x[ii] * F_star_Ux_x + der_coeffs_y[ii] * F_star_Ux_y;
      const auto Uy_ = der_coeffs_x[ii] * F_star_Uy_x + der_coeffs_y[ii] * F_star_Uy_y;
      
      if (! quadrant->is_hanging (ii)){
        increment [ordh  (quadrant->gt (ii))] += h_;
        increment [ordUx (quadrant->gt (ii))] += Ux_;
        increment [ordUy (quadrant->gt (ii))] += Uy_;

      } else {
        
        increment [ordh  (quadrant->gparent(0,ii))] += h_;
        increment [ordh  (quadrant->gparent(1,ii))] += h_;
        
        increment [ordUx (quadrant->gparent(0,ii))] += Ux_;
        increment [ordUx (quadrant->gparent(1,ii))] += Ux_;
        
        increment [ordUy (quadrant->gparent(0,ii))] += Uy_;
        increment [ordUy (quadrant->gparent(1,ii))] += Uy_;
        

      }
    }
    
    
    
    
//    // increment
//    for (ii = 0; ii < 4; ++ii)
//    {
//      loc_incrh [ii] = der_coeffs_x[ii] * F_star_h_x  + der_coeffs_y[ii] * F_star_h_y;
//      loc_incrUx[ii] = der_coeffs_x[ii] * F_star_Ux_x + der_coeffs_y[ii] * F_star_Ux_y;
//      loc_incrUy[ii] = der_coeffs_x[ii] * F_star_Uy_x + der_coeffs_y[ii] * F_star_Uy_y;
//
////      loc_incrh [ii] = der_coeffs_x[ii] * F_star_h_x  + der_coeffs_y[ii] * F_star_h_y;
////      loc_incrUx[ii] = der_coeffs_x[ii] * F_star_Ux_x + der_coeffs_y[ii] * F_star_Ux_y;
////      loc_incrUy[ii] = der_coeffs_x[ii] * F_star_Uy_x + der_coeffs_y[ii] * F_star_Uy_y;
//    }
    

  }
  
  /*
  void
  halfstep_function () {
    
     int jj, kk;
     double tmpdh = 0, tmpdUx = 0, tmpdUy = 0;
     
     loc_midh  = 0.;
     loc_midUx = 0.;
     loc_midUy = 0.;
     
     
     const double hm  = .25 * std::accumulate (hdof.begin(),  hdof.end(),  0.0);
     const double Uxm = .25 * std::accumulate (Uxdof.begin(), Uxdof.end(), 0.0);
     const double Uym = .25 * std::accumulate (Uydof.begin(), Uydof.end(), 0.0);
     
     if (hm > 0.) {
     
     const double dtoptx = Dx / (std::abs (Uxm / hm) + std::sqrt (grav * hm));
     const double dtopty = Dy / (std::abs (Uym / hm) + std::sqrt (grav * hm));
     const double dtopt = dtoptx > dtopty ? dtopty : dtoptx;
     
     if (dt > dtopt) set_dt (dtopt);
     
     for (jj = 0; jj < 4; ++jj) {
     for (kk = 0; kk < 4; ++kk) {
     tmpdh   += wq[kk] * (- shgx[jj][kk] * h_flux_formula_x  (hdof[jj], Uxdof[jj], Uydof[jj])
     - shgy[jj][kk] * h_flux_formula_y  (hdof[jj], Uxdof[jj], Uydof[jj])
     + shp[jj][kk]  * h_src_formula     (hdof[jj], Uxdof[jj], Uydof[jj]));
     tmpdUx  += wq[kk] * (- shgx[jj][kk] * Ux_flux_formula_x (hdof[jj], Uxdof[jj], Uydof[jj])
     - shgy[jj][kk] * Ux_flux_formula_y (hdof[jj], Uxdof[jj], Uydof[jj])
     + shp[jj][kk]  * Ux_src_formula    (dZdx[kk], hdof[jj], Uxdof[jj], Uydof[jj]));
     tmpdUy  += wq[kk] * (- shgx[jj][kk] * Uy_flux_formula_x (hdof[jj], Uxdof[jj], Uydof[jj])
     - shgy[jj][kk] * Uy_flux_formula_y (hdof[jj], Uxdof[jj], Uydof[jj])
     + shp[jj][kk]  * Uy_src_formula    (dZdy[kk], hdof[jj], Uxdof[jj], Uydof[jj]));
     }
     }
     
     loc_midh  = hm  + .5 * dtopt * tmpdh  / area;
     if (loc_midh > 0.) {
     loc_midUx = Uxm + .5 * dtopt * tmpdUx / area;
     loc_midUy = Uym + .5 * dtopt * tmpdUy / area;
     } else
     loc_midh = 0.;
     
     }
    
  }
  
  void
  src_function () {
    loc_srch  = 0.;
    loc_srcUx = 0.;
    loc_srcUy = 0.;
    for (int kk = 0; kk < 4; ++kk) {
      loc_srcUx += wq[kk] * Ux_src_formula (dZdx[kk], loc_midh, loc_midUx, loc_midUy) / area;
      loc_srcUy += wq[kk] * Uy_src_formula (dZdy[kk], loc_midh, loc_midUx, loc_midUy) / area;
    }
  }
  
  
  void
  flux_function () {
    loc_fluxh_x  = h_flux_formula_x  (loc_midh, loc_midUx, loc_midUy);
    loc_fluxh_y  = h_flux_formula_y  (loc_midh, loc_midUx, loc_midUy);
    loc_fluxUx_x = Ux_flux_formula_x (loc_midh, loc_midUx, loc_midUy);
    loc_fluxUx_y = Ux_flux_formula_y (loc_midh, loc_midUx, loc_midUy);
    loc_fluxUy_x = Uy_flux_formula_x (loc_midh, loc_midUx, loc_midUy);
    loc_fluxUy_y = Uy_flux_formula_y (loc_midh, loc_midUx, loc_midUy);
  }
  
  
  void
  incr_function () {
    int ii, jj, kk;
    loc_incrh  = {0, 0, 0, 0};
    loc_incrUx = {0, 0, 0, 0};
    loc_incrUy = {0, 0, 0, 0};
    
    if (loc_midh > 0.)
      for (ii = 0; ii < 4; ++ii)
    for (kk = 0; kk < 4; ++kk) {
      loc_incrh [ii] += wq[kk] * (shp[ii][kk] * loc_srch  + shgx[ii][kk] * loc_fluxh_x  + shgy[ii][kk] * loc_fluxh_y);
      loc_incrUx[ii] += wq[kk] * (shp[ii][kk] * loc_srcUx + shgx[ii][kk] * loc_fluxUx_x + shgy[ii][kk] * loc_fluxUx_y);
      loc_incrUy[ii] += wq[kk] * (shp[ii][kk] * loc_srcUy + shgx[ii][kk] * loc_fluxUy_x + shgy[ii][kk] * loc_fluxUy_y);
    }
    
  }*/
  
  
  void
  set_dt (const double dt_)
  { dt = dt_; }
  
  double
  get_dt ()
  { return dt; }
  
  
  double dt;
  
  double Dx, Dy, area;
  
  
  ///  quadrant vertex (dofs) coordinates
  ///  The assumed numbering for quadrant nodes is
  ///  the following :
  ///
  ///   2------------------3
  ///   |                  |
  ///   |                  |
  ///   |                  |
  ///   |                  |
  ///   0------------------1
  
  std::array<double, 4> xn = {0, 0, 0, 0};
  std::array<double, 4> yn = {0, 0, 0, 0};
  
  // local dofs for state vector components
  std::array<double, 4> hdof  = {0, 0, 0, 0};
  std::array<double, 4> Uxdof = {0, 0, 0, 0};
  std::array<double, 4> Uydof = {0, 0, 0, 0};
  
  
  
  // local buffers for state and flux update
  std::array<double, 4> loc_incrh  = {0, 0, 0, 0};
  std::array<double, 4> loc_incrUx = {0, 0, 0, 0};
  std::array<double, 4> loc_incrUy = {0, 0, 0, 0};
  
  
  
  
  // flux functions
  double
  h_flux_formula_x (const double& h, const double& Ux, const double& Uy)
  { return Ux; }

  double
  h_flux_formula_y (const double& h, const double& Ux, const double& Uy)
  { return Uy; }

  double
  Ux_flux_formula_x (const double& h, const double& Ux, const double& Uy)
  { return h > 0. ?  Ux*Ux/h + grav*h*h/2. : 0.; }

  double
  Ux_flux_formula_y (const double& h, const double& Ux, const double& Uy)
  { return h > 0. ? Uy*Ux/h : 0.; }

  double
  Uy_flux_formula_x (const double& h, const double& Ux, const double& Uy)
  { return h > 0. ? Uy*Ux/h : 0.; }

  double
  Uy_flux_formula_y (const double& h, const double& Ux, const double& Uy)
  { return h > 0. ? Uy*Uy/h + grav*h*h/2. : 0.; }


  // source terms
  double
  h_src_formula (const double& h, const double& Ux, const double& Uy)
  { return (0.); }

  double
  Ux_src_formula (const double& dZdx, const double& h, const double& Ux, const double& Uy)
  { return (-grav*h*dZdx); }

  double
  Uy_src_formula (const double& dZdy, const double& h, const double& Ux, const double& Uy)
  { return (-grav*h*dZdy); }

};








// Re-Define tic and toc to add an MPI_Barrier
#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }


int
main (int argc, char **argv)
{
  // Management of solutions ordering
  ordering ordh  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 0> (gt); };
  ordering ordUx = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 1> (gt); };
  ordering ordUy = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 2> (gt); };

  
  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  
  
  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  TIC ();
  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);
  
  tmsh.set_refine_marker (hanging_refinement);
  tmsh.refine (recursive, 1);
  TOC ("Uniform refinement");
  
  // ln_nodes sono i dof non gli hanging node!! (sono esclusi dal calcolo)
  tmesh::idx_t gn_nodes    = tmsh.num_global_nodes ();
  tmesh::idx_t ln_nodes    = tmsh.num_owned_nodes ();
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();
  
  /// Allocate initial data container
  Q1 sol  (ln_nodes * 3);
  Q1 incr (ln_nodes * 3);
  sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);
  incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);
  
  Q1 mass (ln_nodes * 3);
  bim2a_mass_vector (tmsh, mass, ordh);
  bim2a_mass_vector (tmsh, mass, ordUx);
  bim2a_mass_vector (tmsh, mass, ordUy);
  mass.assemble ();
  // .assemble serve probabilmente in un ambiente parallelo per creare la matrice di massa globale
  
//   gioca un po' con questo per capire la matrice di massa lumpata e come integrare il termine di flusso
//  for (const auto& it : mass.get_owned_data ())
//  {
//    std::cout << it << std::endl;
//  }
//
//  std::cout << "STOP!" << " " << ln_nodes << " " << ln_elements << " " << mass.get_owned_data ().size() <<std::endl;
//  return 0;
  
  
  Q1 Z    (ln_nodes);
  Z.get_owned_data ().assign (Z.get_owned_data ().size (), 0.0);
  

  TIC();
  octave_io_mode m_in = gz_read_mode, m_out = gz_read_mode;
  octave_value v;
  
  octave_io_open (LOADFILENAME_1, m_in, &m_out);
  octave_load (VARNAME_1, v);
  Matrix M = v.matrix_value ();
  dem.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), dem.begin ());
  
  octave_io_open (LOADFILENAME_2, m_in, &m_out);
  octave_load (VARNAME_2, v);
  M = v.matrix_value ();
  basin_mask.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), basin_mask.begin ());
  
  TOC("Load data matrix");
  
  
  // Buffer for export filename
  char filename[255]="";
  
  
  // Initialize initial datas
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
  {
    
    for (int ii = 0; ii < 4; ++ii)
    {
      if (! quadrant->is_hanging (ii)){
        double xx=quadrant->p(0,ii);
        double yy=quadrant->p(1,ii);
        
        sol [ordh (quadrant->gt (ii))] = h0_fun (xx, yy);
        sol [ordUx(quadrant->gt (ii))] = Ux0_fun (xx, yy);
        sol [ordUy(quadrant->gt (ii))] = Uy0_fun (xx, yy);
        
        Z[quadrant->gt (ii)] = dem[global_coord_2_raster(xx,yy)[0]]*0;
        
//        if (h0_fun (xx, yy) !=0) std::cout << h0_fun (xx, yy) << std::endl;
      }
      
      else
      {
        // touch parent nodes to set up distributed vector structure
        
        sol [ordh (quadrant->gparent(0,ii))] += 0.;
        sol [ordh (quadrant->gparent(1,ii))] += 0.;
        sol [ordUx(quadrant->gparent(0,ii))] += 0.;
        sol [ordUx(quadrant->gparent(1,ii))] += 0.;
        sol [ordUy(quadrant->gparent(0,ii))] += 0.;
        sol [ordUy(quadrant->gparent(1,ii))] += 0.;
        
        Z   [quadrant->gparent(0,ii)] += 0.;
        Z   [quadrant->gparent(1,ii)] += 0.;
      }
    }
  }
  // bim2a_solution_with_ghosts in quad_operators.cpp
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy);
  
  bim2a_solution_with_ghosts (tmsh, Z, replace_op);
  
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUy);
  TOC ("compute initial condition");
  
  // create results dir
  if (rank == 0)
  {
    const std::string bashCommand = std::string ( "mkdir -p results" );
    std::system ( bashCommand.c_str() );
  }
  
  
  // Save initial conditions
  sprintf(filename, "results/swe_h_%4.4d", 0);
  tmsh.octbin_export (filename, sol, ordh);
  sprintf(filename, "results/swe_Ux_%4.4d", 0);
  tmsh.octbin_export (filename, sol, ordUx);
  sprintf(filename, "results/swe_Uy_%4.4d", 0);
  tmsh.octbin_export (filename, sol, ordUy);
  sprintf(filename, "results/swe_Z_%4.4d", 0);
  tmsh.octbin_export (filename, Z);
  

  std::vector<double> full_time_vector;
  full_time_vector.reserve (static_cast<int> (T/DELTAT));
  std::vector<double> save_time_vector;
  save_time_vector.reserve (static_cast<int> (T/SAVEDT));
  
  
  
  TG2_rusanov stp(sol, ordh, ordUx, ordUy, Z);

  
  
  // Time loop
  double time = 0.0;
  double deltat = DELTAT;
  if(rank==0) {
    full_time_vector.push_back (0.0);
    save_time_vector.push_back (0.0);
  }
  int count = 0;
  
  double savecount = 0.0;
  
  
  while (time <= T)
  {
    // Reset increment
    TIC();
    incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);
    incr.assemble (replace_op);
    TOC("Reset");
    
    
    TIC();
    stp.set_dt (DELTAT);
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.set_quadrant (quadrant);
      stp.step_function (incr, quadrant);
      
      
//      assemble_vector (quadrant, stp.loc_incrh, incr, ordh);
//      assemble_vector (quadrant, stp.loc_incrUx, incr, ordUx);
//      assemble_vector (quadrant, stp.loc_incrUy, incr, ordUy);
    }

    incr.assemble ();
    TOC("Compute step");
    
    deltat = REDCDT * stp.dt;
    if(rank==0)
      std::cout << "TIME = " << time << ", next dt = " << deltat << std::endl;
    
    MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&deltat), 1, MPI_INT, MPI_MIN, tmsh.comm);
    time += deltat;
    savecount += deltat;
    MPI_Bcast (static_cast<void*> (&time), 1, MPI_DOUBLE, 0, tmsh.comm);
    MPI_Bcast (static_cast<void*> (&savecount), 1, MPI_DOUBLE, 0, tmsh.comm);
    MPI_Barrier (tmsh.comm);
    
    // Print curent time
    if(rank==0)
    {
      std::cout << "TIME = " << time << ", last dt = " << deltat << std::endl;
      full_time_vector.push_back (time);
    }
    
    
    TIC();
    for (auto kk = 0; kk < incr.get_owned_data ().size (); kk++)
    {
//      std::cout << sol.get_owned_data ()[kk] << " ";
      
      
      // this is the increment, sol.get_owned_data () is a vector probably because it has .assign function
      sol.get_owned_data ()[kk] += deltat * incr.get_owned_data ()[kk] / mass.get_owned_data ()[kk];
      
      
      
//      std::cout << incr.get_owned_data ()[kk] << " " << sol.get_owned_data ()[kk] << std::endl;
    }
//    return 0;

    
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant) {
      
      for (int i = 0; i < 4; ++i) {
        if (! quadrant->is_hanging (i)) {
          
          auto boundary_idx = quadrant->e (i);
          auto boundary_idxx = quadrant->ex (i);
          auto boundary_idxy = quadrant->ey (i);
          
          if (boundary_idxx != tmesh::quadrant_t::NOT_ON_BOUNDARY) {
            sol[ordUx (quadrant->gt (i))] = 0;
//            sol[ordUy (quadrant->gt (i))] = 0;
            
//            std::cout << boundary_idx << std::endl;
            
          }
          
          if (boundary_idxy != tmesh::quadrant_t::NOT_ON_BOUNDARY) {
            sol[ordUy (quadrant->gt (i))] = 0;
              //            sol[ordUy (quadrant->gt (i))] = 0;
            
              //            std::cout << boundary_idx << std::endl;
            
          }
          
        }
      }
    }
//    return 0;
    /*
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant) {
      auto tree_idx = quadrant->get_tree_idx ();
      
      for (int i = 0; i < 4; ++i) {
        if (! quadrant->is_hanging (i)) {
          
          auto boundary_idx = quadrant->e (i);
          auto boundary_idxx = quadrant->ex (i);
          auto boundary_idxy = quadrant->ey (i);
          
          // If current node is on boundary
          if (boundary_idx != tmesh::quadrant_t::NOT_ON_BOUNDARY) {
            
            // Loop over all the boundary conditions on h.
            for (size_t bc = 0; bc < bcsh.size (); ++bc)
              // If this boundary condition matches with
              // the current node.
            if (std::get<0> (bcsh[bc]) == tree_idx
                && (std::get<1> (bcsh[bc]) == boundary_idx
                    || std::get<1> (bcsh[bc]) == boundary_idxx
                    || std::get<1> (bcsh[bc]) == boundary_idxy)) {
                // Evaluate bc at current node.
              sol[ordh (quadrant->gt (i))] = (std::get<2> (bcsh[bc])) (quadrant->p (0, i), quadrant->p (1, i));
            }
            
            // Loop over all the boundary conditions on Ux.
            for (size_t bc = 0; bc < bcsUx.size (); ++bc)
              // If this boundary condition matches with
              // the current node.
            if (std::get<0> (bcsUx[bc]) == tree_idx
                && (std::get<1> (bcsUx[bc]) == boundary_idx
                    || std::get<1> (bcsUx[bc]) == boundary_idxx
                    || std::get<1> (bcsUx[bc]) == boundary_idxy)) {
                // Evaluate bc at current node.
              sol[ordUx (quadrant->gt (i))] = (std::get<2> (bcsUx[bc])) (quadrant->p (0, i), quadrant->p (1, i));
            }
            
            // Loop over all the boundary conditions on Uy.
            for (size_t bc = 0; bc < bcsUy.size (); ++bc)
              // If this boundary condition matches with
              // the current node.
            if (std::get<0> (bcsUy[bc]) == tree_idx
                && (std::get<1> (bcsUy[bc]) == boundary_idx
                    || std::get<1> (bcsUy[bc]) == boundary_idxx
                    || std::get<1> (bcsUy[bc]) == boundary_idxy)) {
                // Evaluate bc at current node.
              sol[ordUy (quadrant->gt (i))] = (std::get<2> (bcsUy[bc])) (quadrant->p (0, i), quadrant->p (1, i));
            }
          }
        }
      }
    }*/
    sol.assemble (replace_op);
    TOC("Apply increment");
    
      // Save solution
    if (savecount >= SAVEDT) {
      TIC();
      if (rank == 0)
        std::cout << "savecount = " << savecount << std::endl;
      count++;
      save_time_vector.push_back (time);
      sprintf(filename, "results/swe_h_%4.4d",   count);
      tmsh.octbin_export (filename, sol, ordh);
      sprintf(filename, "results/swe_Ux_%4.4d",  count);
      tmsh.octbin_export (filename, sol, ordUx);
      sprintf(filename, "results/swe_Uy_%4.4d",  count);
      tmsh.octbin_export (filename, sol, ordUy);
      sprintf(filename, "results/swe_Z_%4.4d",  count);
      tmsh.octbin_export (filename, Z);
      savecount = 0.0;
      TOC("Exporting solution");
    }
    
    
  }
  
  
  if (rank == 0)
  {
    ColumnVector vtmp (save_time_vector.size ());
    std::copy (save_time_vector.begin (), save_time_vector.end (), vtmp.fortran_vec ());
    octave_io_mode m;
    octave_io_open ("timesteps.octbin", gz_write_mode, &m);
    octave_save ("save_time", vtmp);
    vtmp.resize (full_time_vector.size ());
    std::copy (full_time_vector.begin (), full_time_vector.end (), vtmp.fortran_vec ());
    octave_save ("full_time", vtmp);
    octave_io_close ();
  }
  
  
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }
  
  MPI_Finalize ();
  
  
  return 0;
  
}
