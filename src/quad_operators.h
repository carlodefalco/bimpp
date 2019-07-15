#ifndef HAVE_QUAD_OPERATORS_H
#define HAVE_QUAD_OPERATORS_H 1

#include <tmesh.h>
#include <bim_distributed_vector.h>
#include <bim_ordering.h>
#include <bim_sparse.h>
#include <operators.h>

#include <functional>
#include <tuple>
#include <vector>

/// f(x, y).
using func = std::function<double (double, double)>;

/// f(quadrant, node index).
using func_quad = std::function<double (tmesh::quadrant_iterator,
                                        tmesh::idx_t)>;

/// Tree index, boundary index, function.
using dirichlet_bcs = std::vector<std::tuple<int, int, func>>;
using dirichlet_bcs_quad = std::vector<std::tuple<int, int, func_quad>>;

template <class T>
using q1_vec = T; // std::vector<double> or distributed_vector.

template <class T>
using gradient = std::pair<T, T>;

/// Nodes, faces, cell midpoint dofs.
using q2_vec = std::vector<std::array<double, 9>>;

/// Function to mark if a quadrant has to be taken into
/// account when computing the recovered gradient.
using active_fun = std::function<bool (tmesh::quadrant_iterator)>;

void
bim2a_structure (tmesh &tmsh,
                 sparse_matrix& A,
                 const ordering& ordr = default_ord,
                 const ordering& ordc = default_ord);

void
bim2a_laplacian (tmesh & mesh,
                 const std::vector<double>& alpha,
                 sparse_matrix& A,
                 const ordering& ordr = default_ord,
                 const ordering& ordc = default_ord);

template <class T>
void
bim2a_laplacian_eafe (tmesh& mesh,
                      const std::vector<double>& D,
                      const T& alpha,
                      sparse_matrix& A,
                      const ordering& ordr,
                      const ordering& ordc);

template <class T>
void
bim2a_advection_diffusion (tmesh & mesh,
                           const std::vector<double>& alpha,
                           const T& psi,
                           sparse_matrix& A,
                           bool symmetric = true,
                           const ordering& ordr = default_ord,
                           const ordering& ordc = default_ord);

template <class T>
void
bim2a_advection_upwind (tmesh & mesh,
                        const T& psi,
                        sparse_matrix& A,
                        const ordering& ordr = default_ord,
                        const ordering& ordc = default_ord);

template <class T>
void
bim2a_advection_eafe_diffusion (tmesh & mesh,
                                const T& alpha,
                                const T& psi,
                                sparse_matrix& A,
                                const ordering& ordr = default_ord,
                                const ordering& ordc = default_ord);

template <class T>
void
bim2a_reaction (tmesh& mesh,
                const std::vector<double>& delta,
                const T& zeta,
                sparse_matrix& A,
                const ordering& ordr = default_ord,
                const ordering& ordc = default_ord);

template <class T>
void
bim2a_rhs (tmesh& mesh,
           const std::vector<double>& f,
           const T& g,
           T& rhs,
           const ordering& ord = default_ord);

void
bim2a_solution_with_ghosts (tmesh& mesh,
                            distributed_vector& rhs,
                            const binary_operator &op = std::plus<double> (),
                            const ordering& ord = default_ord,
                            bool remap_assemble = true);

template <class T>
void
bim2a_boundary_mass (tmesh & mesh,
                     const int & tree_idx,
                     const int & boundary_idx,
                     T & M,
                     const func_quad & fun =
                     [] (tmesh::quadrant_iterator, tmesh::idx_t)
                       {return 1;},
                     const ordering & ord = default_ord);

/// rhs_or_A:
/// 0, for both matrix and rhs;
/// 1, for only rhs;
/// 2, for only matrix.
template <class T>
void
bim2a_dirichlet_bc (tmesh& mesh, const dirichlet_bcs& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord = default_ord,
                    const bool& only_rhs = false);

template <class T>
void
bim2a_dirichlet_bc (tmesh& mesh, const dirichlet_bcs_quad& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord = default_ord,
                    const bool& only_rhs = false);

template <class T>
void
bim2a_robin_bc_loc (tmesh& mesh,
                    sparse_matrix& A,
                    T& rhs,
                    const double& M_boundary_loc,
                    const unsigned int& row,
                    const double& value_A,
                    const double& value_rhs,
                    const bool& only_rhs = false,
		    int start = -1,
		    int end = -1);

template <class T>
void
interpolate_vector (tmesh & mesh,
                    T & vec_in,
                    T & vec_out,
                    const ordering & ord = default_ord);





template <class T>
gradient<T>
bim2c_quadtree_pde_recovered_gradient (tmesh & mesh,
                                       const T & u,
                                       active_fun is_active =
                                       [] (tmesh::quadrant_iterator)
                                         {return true;});

template <class T>
q2_vec
bim2c_quadtree_pde_recovered_solution (tmesh & mesh,
                                       const T & u,
                                       const gradient<T> & du);

/// Compute ||grad^* u - grad u||_L^2(q).
template <class T>
double
estimator_grad (tmesh::quadrant_iterator q,
                const gradient<T> & du_star,
                const T & u);

template <class T>
int
zz_marker_grad (tmesh::quadrant_iterator q,
                const gradient<T> & du_star,
                const T & u,
                double limit);

/// Compute ||u^* - u||_L^2(q).
template <class T>
double
estimator_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const T & u);

template <class T>
int
zz_marker_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const T & u,
               double limit);

/// Compute ||u - u_ex||_L^2(q).
template <class T>
double
l2_error (tmesh::quadrant_iterator q,
          const func & u_ex,
          const T & u);

/// Compute |u - u_ex|_H^1(q).
template <class T>
double
semih1_error (tmesh::quadrant_iterator q,
              const func & dudx_ex,
              const func & dudy_ex,
              const T & u);

// Compute ||u_star - u_ex||_L^2(q).
double
l2_star_error (tmesh::quadrant_iterator q,
               const func & u_ex,
               const q2_vec & ustar);

/// Compute |du_star - grad(u_ex)|_H^1(q).
template <class T>
double
semih1_star_error (tmesh::quadrant_iterator q,
                   const func & dudx_ex,
                   const func & dudy_ex,
                   const gradient<T> & du_star);


#endif
