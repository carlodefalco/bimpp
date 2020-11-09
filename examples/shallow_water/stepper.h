#include <bim_distributed_vector.h>
#include <tmesh.h>
#include <quad_operators.h>

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
    static_cast<T*>(this)->update_function ();
  }

  void
  set_quadrant (tmesh::quadrant_iterator q)
  { quadrant = q; }

  const Q1& state_vector; 
  const ordering& ordh;
  const ordering& ordUx;
  const ordering& ordUy;

  std::array<double, 4> xn = {0, 0, 0, 0};
  std::array<double, 4> yn[4] = {0, 0, 0, 0};
  
  std::array<double, 4> hdof  = {0, 0, 0, 0};
  std::array<double, 4> Uxdof = {0, 0, 0, 0};
  std::array<double, 4> Uydof = {0, 0, 0, 0};

  std::array<double, 4*nquad> shp;
  std::array<double, 4*nquad> shgx;
  std::array<double, 4*nquad> shgy;

  double xq[nquad];
  double yq[nquad];

  double wq[nquad];
  
  std::array<double, 4> loc_incrh  = {0, 0, 0, 0};
  std::array<double, 4> loc_incrUx = {0, 0, 0, 0};
  std::array<double, 4> loc_incrUy = {0, 0, 0, 0};
  
private :

  tmesh::quadrant_iterator quadrant;

};


class
peraire_stepper
  : public abstract_stepper<peraire_stepper, 4>
{
public :
  using basetype = abstract_stepper<peraire_stepper, 4>;
  peraire_stepper () = delete;
  peraire_stepper  (tmesh& msh,
                    const basetype::Q1& state,
                    const ordering& oh,
                    const ordering& oUx,
                    const ordering& oUy)
    : basetype(msh, state, oh, oUx, oUy)
  { };


  void
  update_function () {
    for (int ii = 0; ii < 4; ++ii){
      loc_incrh [ii] = - 0.01 * hdof [ii];
      loc_incrUx[ii] = - 0.01 * Uxdof[ii];
      loc_incrUy[ii] = - 0.01 * Uydof[ii];
    }
  };
  
};
