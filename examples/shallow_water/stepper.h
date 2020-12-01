#include <numeric>
#include <bim_distributed_vector.h>
#include <tmesh.h>
#include <quad_operators.h>


/// Base interface (template) class
template<class T, int nquad>
class abstract_stepper
{

  
public :

  using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
  using Q0  = std::vector<double>;         // Typedef for local q_0 vector

  abstract_stepper () = delete;
  abstract_stepper (tmesh& msh,
                    const Q1& state,
                    const ordering& oh,
                    const ordering& oUx,
                    const ordering& oUy)
    : quadrant(msh.begin_quadrant_sweep ()), state_vector(state), ordh(oh), ordUx(oUx), ordUy(oUy)
  {  };

  void
  update_halfstep ()
  { static_cast<T*>(this)->halfstep_function (); }

  void
  update_flux ()
  { static_cast<T*>(this)->flux_function (); }

  void
  update_src ()
  { static_cast<T*>(this)->src_function (); }

  void
  update_state_incr ()
  { static_cast<T*>(this)->incr_function (); }

  void
  set_quadrant (tmesh::quadrant_iterator q)
  {
    quadrant = q;
    for (int ii = 0; ii < 4; ++ii) {
      xn[ii] = quadrant->p(0, ii);
      yn[ii] = quadrant->p(1, ii);
    }

    for (int ii = 0; ii < 4; ++ii){
      if (! quadrant->is_hanging (ii)){  
        hdof[ii] = state_vector [ordh (quadrant->gt (ii))];
        Uxdof[ii] = state_vector [ordUx (quadrant->gt (ii))];
        Uydof[ii] = state_vector [ordUy (quadrant->gt (ii))];
      } else {
        hdof[ii]  = .5 * (state_vector [ordh (quadrant->gparent(0,ii))] +
                          state_vector [ordh (quadrant->gparent(1,ii))]);
        Uxdof[ii] = .5 * (state_vector [ordUx (quadrant->gparent(0,ii))] +
                          state_vector [ordUx (quadrant->gparent(1,ii))]);
        Uydof[ii] = .5 * (state_vector [ordUy (quadrant->gparent(0,ii))] +
                          state_vector [ordUy (quadrant->gparent(1,ii))]);
      }
    }
    static_cast<T*>(this)->basis_functions ();
    static_cast<T*>(this)->update_coeffs ();
  }

  void
  set_dt (const double dt_)
  { dt = dt_; }

  double
  get_dt ()
  { return dt; }

  void
  set_flux (const double Fh_x, const double Fh_y, const double FUx_x,
            const double FUx_y, const double FUy_x, const double FUy_y)
  { loc_fluxh_x = Fh_x; loc_fluxUx_x = FUx_x; loc_fluxUy_x = FUy_x;
    loc_fluxh_y = Fh_y; loc_fluxUx_y = FUx_y; loc_fluxUy_y = FUy_y;}

  void
  get_flux (double& Fh_x, double& Fh_y, double& FUx_x,
            double& FUx_y, double& FUy_x, double& FUy_y)
  { Fh_x = loc_fluxh_x; FUx_x = loc_fluxUx_x ; FUy_x = loc_fluxUy_x;
    Fh_y = loc_fluxh_y; FUx_y = loc_fluxUx_y ; FUy_y = loc_fluxUy_y;}

  void
  set_src (const double Sh, const double SUx, const double SUy)
  { loc_srch = Sh; loc_srcUx = SUx; loc_srcUy = SUy; }

  void
  get_src (double& Sh, double& SUx, double& SUy)
  { Sh = loc_srch; SUx = loc_srcUx ; SUy = loc_srcUy; }
  
  const Q1& state_vector; 
  const ordering& ordh;
  const ordering& ordUx;
  const ordering& ordUy;

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

  std::array<std::array<double, nquad>, 4> shp;
  std::array<std::array<double, nquad>, 4> shgx;
  std::array<std::array<double, nquad>, 4> shgy;

  // quadrature nodes
  int get_nquad() { return nquad; }
  std::array<double, nquad> xq;
  std::array<double, nquad> yq;

  // quadrature weights
  std::array<double, nquad> wq;

  // local buffers for state and flux update
  std::array<double, 4> loc_incrh  = {0, 0, 0, 0};
  std::array<double, 4> loc_incrUx = {0, 0, 0, 0};
  std::array<double, 4> loc_incrUy = {0, 0, 0, 0};

  double loc_midh  = 0;
  double loc_midUx = 0;
  double loc_midUy = 0;
  
  double loc_fluxh_x  = 0;
  double loc_fluxh_y  = 0;
  double loc_fluxUx_x = 0;
  double loc_fluxUx_y = 0;
  double loc_fluxUy_x = 0;
  double loc_fluxUy_y = 0;

  double loc_srch  = 0;
  double loc_srcUx = 0;
  double loc_srcUy = 0;

  double dt = 1.0;
  
protected :

  tmesh::quadrant_iterator quadrant;

};


/// Default implementation specializing to a particular quadrature
template<class T>
class
default_stepper_4
  : public abstract_stepper<T, 4>
{

public:
  using  basetype = abstract_stepper<T, 4>;


  T& top() { return static_cast<T&>(*this); }

  default_stepper_4 (tmesh& msh,
                     const typename basetype::Q1& state,
                     const ordering& oh,
                     const ordering& oUx,
                     const ordering& oUy)
    : basetype(msh, state, oh, oUx, oUy)
  { }
  
  void
  basis_functions () {
    compute_quad_points ();
    compute_shp ();
    compute_shgx ();
    compute_shgy ();
  }

  void
  update_coeffs () { }
  
  void
  compute_quad_points () {
    top().xq = top().xn;
    top().yq = top().yn;
    const double hxhyby4 = .25 * (top().xn[1] - top().xn[0]) *
      (top().yn[2] - top().yn[0]);
    top().wq = {hxhyby4, hxhyby4, hxhyby4, hxhyby4};
  }

  void
  compute_shp () {
    top().shp[0] = {1, 0, 0, 0};
    top().shp[1] = {0, 1, 0, 0};
    top().shp[2] = {0, 0, 1, 0};
    top().shp[3] = {0, 0, 0, 1};
  }

  void
  compute_shgx () {
    const double hx = top().xn[1]-top().xn[0]; 
    top().shgx[0] = {-1./hx, -1./hx,  0,      0};
    top().shgx[1] = {1./hx,   1./hx,  0,      0};
    top().shgx[2] = {0,       0,     -1./hx, -1./hx};
    top().shgx[3] = {0,       0,      1./hx,  1./hx};
  }

  void
  compute_shgy () {
    const double hy = top().yn[2]-top().yn[0]; 
    top().shgy[0] = {-1./hy,  0,     -1./hy,  0};
    top().shgy[1] = { 0,     -1./hy,  0,     -1./hy};
    top().shgy[2] = { 1./hy,  0,      1./hy,  0};
    top().shgy[3] = { 0,      1./hy,  0,      1./hy};
  }
  
};



/// Final concrete class implementing a particular model/method
class
peraire_stepper
  : public default_stepper_4<peraire_stepper>
{

private :
  
  const basetype::Q1& Z;
  std::array<double, 4> dZdx  = {0, 0, 0, 0};
  std::array<double, 4> dZdy  = {0, 0, 0, 0};
  
public :
  
  using basetype = default_stepper_4<peraire_stepper>;
  peraire_stepper () = delete;
  peraire_stepper  (tmesh& msh,
                    const basetype::Q1& state,
                    const ordering& oh,
                    const ordering& oUx,
                    const ordering& oUy,
                    const basetype::Q1& elevation)
    : basetype(msh, state, oh, oUx, oUy), Z(elevation)
  { }

  void
  update_coeffs () {

    dZdx.fill (0.0);
    dZdy.fill (0.0);
    for (int kk = 0; kk < get_nquad (); ++kk) {
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
    }

  }
  
  static constexpr double grav = 9.81;
  double
  h_flux_formula_x (double h, double Ux, double Uy)
  { return h > 0. ? Ux : 0.; }

  double
  h_flux_formula_y (double h, double Ux, double Uy)
  { return h > 0. ? Uy : 0.; }

  double
  Ux_flux_formula_x (double h, double Ux, double Uy)
  { return h > 0. ?  Ux*Ux/h + grav*h*h/2. : 0.; }
  
  double
  Ux_flux_formula_y (double h, double Ux, double Uy)
  { return h > 0. ? Uy*Ux/h : 0.; }

  double
  Uy_flux_formula_x (double h, double Ux, double Uy)
  { return h > 0. ? Uy*Ux/h : 0.; }
  
  double
  Uy_flux_formula_y (double h, double Ux, double Uy)
  { return h > 0. ? Uy*Uy/h + grav*h*h/2. : 0.; }

  double
  h_src_formula (double h, double Ux, double Uy)
  { return (0.); }

  double
  Ux_src_formula (double dZdx, double h, double Ux, double Uy)
  { return (-grav*h*dZdx); }

  double
  Uy_src_formula (double dZdy, double h, double Ux, double Uy)
  { return (-grav*h*dZdy); }

   void
  halfstep_function () {
    int jj, kk;
    double tmpdh = 0, tmpdUx = 0, tmpdUy = 0,
      tmph = 0, tmpUx = 0, tmpUy = 0;
    
    double area = (xn[1] - xn[0]) * (yn[2] - yn[0]);

    loc_midh  = 0.;
    loc_midUx = 0.;
    loc_midUy = 0.;

    const double hx = xn[1]-xn[0];
    const double hy = yn[2]-yn[0];
    const double hm  = .25 * std::accumulate (hdof.begin(),  hdof.end(),  0.0);
    const double Uxm = .25 * std::accumulate (Uxdof.begin(), Uxdof.end(), 0.0);
    const double Uym = .25 * std::accumulate (Uydof.begin(), Uydof.end(), 0.0);
      
    if (hm > 0.) {

      const double dtoptx = hx / (std::abs (Uxm / hm) + std::sqrt (grav * hm));
      const double dtopty = hy / (std::abs (Uym / hm) + std::sqrt (grav * hm));
      const double dtopt = dtoptx > dtopty ? dtopty : dtoptx;
    
      if (get_dt () > dtopt) set_dt (dtopt);
      
      for (jj = 0; jj < 4; ++jj) {
        for (kk = 0; kk < get_nquad (); ++kk) {
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
  flux_function () {
      loc_fluxh_x  = h_flux_formula_x  (loc_midh, loc_midUx, loc_midUy);
      loc_fluxh_y  = h_flux_formula_y  (loc_midh, loc_midUx, loc_midUy);
      loc_fluxUx_x = Ux_flux_formula_x (loc_midh, loc_midUx, loc_midUy);
      loc_fluxUx_y = Ux_flux_formula_y (loc_midh, loc_midUx, loc_midUy);
      loc_fluxUy_x = Uy_flux_formula_x (loc_midh, loc_midUx, loc_midUy);
      loc_fluxUy_y = Uy_flux_formula_y (loc_midh, loc_midUx, loc_midUy);
  }

  void
  src_function () {
      loc_srch  = 0.;
      loc_srcUx = 0.;
      loc_srcUy = 0.;
      double area = (xn[1] - xn[0]) * (yn[2] - yn[0]);
      for (int kk = 0; kk < 4; ++kk) {
        loc_srcUx += wq[kk] * Ux_src_formula (dZdx[kk], loc_midh, loc_midUx, loc_midUy) / area;
        loc_srcUy += wq[kk] * Uy_src_formula (dZdy[kk], loc_midh, loc_midUx, loc_midUy) / area;
    }
  }
  
  void
  incr_function () {
    int ii, jj, kk;
    loc_incrh  = {0, 0, 0, 0};
    loc_incrUx = {0, 0, 0, 0};
    loc_incrUy = {0, 0, 0, 0};

    if (loc_midh > 0.) 
      for (ii = 0; ii < 4; ++ii) 
          for (kk = 0; kk < get_nquad (); ++kk) {
            loc_incrh [ii] += wq[kk] * (shp[ii][kk] * loc_srch  + shgx[ii][kk] * loc_fluxh_x  + shgy[ii][kk] * loc_fluxh_y);
            loc_incrUx[ii] += wq[kk] * (shp[ii][kk] * loc_srcUx + shgx[ii][kk] * loc_fluxUx_x + shgy[ii][kk] * loc_fluxUx_y);
            loc_incrUy[ii] += wq[kk] * (shp[ii][kk] * loc_srcUy + shgx[ii][kk] * loc_fluxUy_x + shgy[ii][kk] * loc_fluxUy_y);
          }
          
  }


};


