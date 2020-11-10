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
  update ()
  {
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

    for (int ii = 0; ii < 4; ++ii) {
      xn[ii] = quadrant->p(0, ii);
      yn[ii] = quadrant->p(1, ii);
    }

    static_cast<T*>(this)->basis_functions ();
    static_cast<T*>(this)->update_function ();
    
  }

  void
  set_quadrant (tmesh::quadrant_iterator q)
  { quadrant = q; }

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
  
  std::array<double, 4> loc_incrh  = {0, 0, 0, 0};
  std::array<double, 4> loc_incrUx = {0, 0, 0, 0};
  std::array<double, 4> loc_incrUy = {0, 0, 0, 0};
  
private :

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
  compute_quad_points () {
    top().xq = top().xn;
    top().yq = top().yn;
    double hxhyby4 = .25 * (top().xn[1] - top().xn[0]) * (top().yn[2] - top().yn[0]);
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
  
public :
  
  using basetype = default_stepper_4<peraire_stepper>;
  peraire_stepper () = delete;
  peraire_stepper  (tmesh& msh,
                    const basetype::Q1& state,
                    const ordering& oh,
                    const ordering& oUx,
                    const ordering& oUy)
    : basetype(msh, state, oh, oUx, oUy)
  { }

  void
  update_function () {
    int ii, jj, kk;
    loc_incrh  = {0, 0, 0, 0};
    loc_incrUx = {0, 0, 0, 0};
    loc_incrUy = {0, 0, 0, 0};
    
    for (ii = 0; ii < 4; ++ii) {
      for (jj = 0; jj < 4; ++jj) {
        for (kk = 0; kk < get_nquad (); ++kk) {
          loc_incrh [ii] += - 2.e1 * wq[kk] * (shgx[ii][kk] * shgx[jj][kk] +
                                               shgy[ii][kk] * shgy[jj][kk]) * hdof[jj];
          loc_incrUx[ii] += - 2.e1 * wq[kk] * (shgx[ii][kk] * shgx[jj][kk] +
                                               shgy[ii][kk] * shgy[jj][kk]) * Uxdof[jj];
          loc_incrUy[ii] += - 2.e1 * wq[kk] * (shgx[ii][kk] * shgx[jj][kk] +
                                               shgy[ii][kk] * shgy[jj][kk]) * Uydof[jj];
        }
      }
    }
  }


};
