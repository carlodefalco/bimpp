#ifndef HAVE_QUAD_OPERATORS_3D_H
#define HAVE_QUAD_OPERATORS_3D_H 1

#include "bimutils.h"
#include <tmesh_3d.h>
#include <bim_distributed_vector.h>
#include <bim_ordering.h>
#include <bim_sparse.h>

#include <cmath>
#include <functional>
#include <tuple>
#include <vector>

// CDF : FIXME : ORDERINGS ARE PASSED BUT NOT ACTUALLY USED !!!

/// f(x, y, z).
using func3 = std::function<double (double, double, double)>;

/// f(quadrant, node index).
using func3_quad = std::function<double (tmesh_3d::quadrant_iterator,
					 tmesh_3d::idx_t)>;

/// Tree index, boundary index, function.
using dirichlet_bcs3 = std::vector<std::tuple<int, int, func3>>;
using dirichlet_bcs3_quad = std::vector<std::tuple<int, int, func3_quad>>;

template <class T>
using q1_vec = T;   // std::vector<double> or distributed_vector

template <class T>
using gradient3 = std::tuple<T, T, T>;

/// Nodes, faces, cell midpoint dofs.
using q2_vec3 = std::vector<std::array<double, 27>>;

/// Function to mark if a quadrant has to be taken into
/// account when computing the recovered gradient.
using active_fun3 = std::function<bool (tmesh_3d::quadrant_iterator)>;


void
bim3a_structure (tmesh_3d &tmsh,
                 sparse_matrix& A,
                 const ordering& ordr = default_ord,
                 const ordering& ordc = default_ord);

void
bim3a_laplacian (tmesh_3d & mesh,
                 const std::vector<double>& alpha,
                 sparse_matrix& A,
                 const ordering& ordr = default_ord,
                 const ordering& ordc = default_ord);

// template <class T>
// void
// bim3a_laplacian_eafe (tmesh_3d& mesh,
//                       const std::vector<double>& D,
//                       const T& alpha,
//                       sparse_matrix& A,
//                       const ordering& ordr,
//                       const ordering& ordc);

void
bim3a_laplacian_eafe (tmesh_3d & mesh,
                      distributed_vector& alpha,
                      sparse_matrix& A,
                      const ordering& ordr = default_ord,
                      const ordering& ordc = default_ord);

// void
// bim3a_laplacian_frac (tmesh_3d & mesh,
//                       distributed_vector& alpha,
//                       sparse_matrix& A,
//                       std::function<std::array<double,12> (double, double, double, double, double, double)> fract,
//                       const ordering& ordr = default_ord,
//                       const ordering& ordc = default_ord);

void
bim3a_laplacian_frac (tmesh_3d & mesh,
                      distributed_vector& alpha,
                      sparse_matrix& A,
                      std::function<std::array<double,12> (tmesh_3d::quadrant_iterator&)> fract,
                      const ordering& ordr = default_ord,
                      const ordering& ordc = default_ord);

template <class T>
void
bim3a_advection_diffusion (tmesh_3d & mesh,
                           const std::vector<double>& alpha,
                           const T& psi,
                           sparse_matrix& A,
                           bool symmetric = true,
                           const ordering& ordr = default_ord,
                           const ordering& ordc = default_ord);

template <class T>
void
bim3a_advection_eafe_diffusion (tmesh_3d & mesh,
                                const T& alpha,
                                const T& psi,
                                sparse_matrix& A,
                                const ordering& ordr = default_ord,
                                const ordering& ordc = default_ord);

template <class T>
void
bim3a_reaction (tmesh_3d& mesh,
                const std::vector<double>& delta,
                const T& zeta,
                sparse_matrix& A,
                const ordering& ordr = default_ord,
                const ordering& ordc = default_ord);


void
bim3a_reaction_frac (tmesh_3d& mesh,
                const distributed_vector& delta,
                const distributed_vector& zeta,
                sparse_matrix& A,
                std::function<std::array<double,12> (tmesh_3d::quadrant_iterator&)> fract,
                const ordering& ordr = default_ord,
                const ordering& ordc = default_ord);

template <class T>
void
bim3a_rhs (tmesh_3d& mesh,
           const std::vector<double>& f,
           const T& g, T& rhs,
           const ordering& ord = default_ord);

void
bim3a_rhs_frac (tmesh_3d& mesh,
                const distributed_vector& f,
                const distributed_vector& g, 
                distributed_vector& rhs,
                std::function<std::array<double,12> (tmesh_3d::quadrant_iterator&)> fract,
                const ordering& ord= default_ord);

void
bim3a_solution_with_ghosts (tmesh_3d& mesh,
                            distributed_vector& rhs,
                            const binary_operator &op = std::plus<double> (),
                            const ordering& ord = default_ord,
                            bool remap_assemble = true);

template <class T>
void
bim3a_boundary_mass (tmesh_3d & mesh,
                     const int & tree_idx,
                     const int & boundary_idx,
                     T & M,
                     const func3_quad & fun =
                     [] (tmesh_3d::quadrant_iterator, tmesh_3d::idx_t)
                       {return 1;},
                     const ordering & ord = default_ord);

/// rhs_or_A:
/// 0, for both matrix and rhs;
/// 1, for only rhs;
/// 2, for only matrix.
template <class T>
void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord = default_ord,
                    const bool& only_rhs = false);

template <class T>
void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3_quad& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord = default_ord,
                    const bool& only_rhs = false);

template <class T>
void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ordr,
                    const ordering& ordc,
                    const bool& only_rhs = false);

// template <class T>
// void
// interpolate_vector (tmesh_3d & mesh,
//                     T & vec_in,
//                     T & vec_out,
//                     const ordering & ord = default_ord);



 double
dudx (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u);
 double
dudy (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u);
 double
dudz (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u);













template <class T>
static double
nedelec_gradient (tmesh_3d::quadrant_iterator & q,
                  const T& u, size_t i);

template <class T>
gradient3<T>
bim3c_quadtree_pde_recovered_gradient (tmesh_3d& mesh,
                                       const T& u,
                                       active_fun3 is_active =
                                        [] (tmesh_3d::quadrant_iterator)
                                         {return true;});

template <class T>
q2_vec3
bim3c_quadtree_pde_recovered_solution (tmesh_3d& mesh,
                                       const T& u,
                                       const gradient3<T>& du);

/// Compute ||grad^* u - grad u||_L^2(q).
template <class T>
double
estimator_grad (tmesh_3d::quadrant_iterator q,
                const gradient3<T>& du_star,
                const T& u);

template <class T>
int
zz_marker_grad (tmesh_3d::quadrant_iterator q,
                const gradient3<T>& du_star,
                const T& u,
                double limit);

/// Compute ||u^* - u||_L^2(q).
template <class T>
double
estimator_sol (tmesh_3d::quadrant_iterator q,
               const q2_vec3 & ustar,
               const T & u);

template <class T>
int
zz_marker_sol (tmesh_3d::quadrant_iterator q,
               const q2_vec3& ustar,
               const T& u,
               double limit);

/// Compute ||u - u_ex||_L^2(q).
template <class T>
double
l2_error (tmesh_3d::quadrant_iterator q,
          const func3 & u_ex,
          const T & u);

/// Compute |u - u_ex|_H^1(q).
template <class T>
double
semih1_error (tmesh_3d::quadrant_iterator q,
              const func3 & dudx_ex,
              const func3 & dudy_ex,
              const func3 & dudz_ex,
              const T & u);

/// Compute ||u_star - u_ex||_L^2(q).
double
l2_star_error (tmesh_3d::quadrant_iterator q,
               const func3 & u_ex,
               const q2_vec3 & ustar);

/// Compute ||du_star - grad(u_ex)||_L^2(q).
template <class T>
double
semih1_star_error (tmesh_3d::quadrant_iterator q,
                   const func3 & dudx_ex,
                   const func3 & dudy_ex,
                   const func3 & dudz_ex,
                   const gradient3<T> & du_star);

#endif
