#ifndef HAVE_QUAD_OPERATORS_H
#define HAVE_QUAD_OPERATORS_H 1

#include "bim_sparse.h"
#include "operators.h"
#include "tmesh.h"

#include <functional>
#include <tuple>
#include <vector>

/// f(x, y).
using func = std::function<double (double, double)>;

/// f(quadrant, node index).
using func_quad = std::function<double (tmesh::quadrant_iterator, tmesh::idx_t)>;

/// Tree index, boundary index, function.
using dirichlet_bcs = std::vector<std::tuple<int, int, func>>;
using dirichlet_bcs_quad = std::vector<std::tuple<int, int, func_quad>>;

using q1_vec = std::vector<double>;

using gradient = std::pair<q1_vec, q1_vec>;

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

template <class T>
void
bim2a_advection_diffusion (tmesh & mesh,
                           const std::vector<double>& alpha,
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

template <class T>
void
bim2a_boundary_mass (tmesh & mesh,
		     const int & tree_idx,
		     const int & boundary_idx,
		     T & M,
		     const func_quad & fun =
		     [] (tmesh::quadrant_iterator, tmesh::idx_t)
		       {return 1;}
		     );

template <class T>
void
bim2a_dirichlet_bc (tmesh& mesh, const dirichlet_bcs& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord = default_ord);

template <class T>
void
bim2a_dirichlet_bc (tmesh& mesh, const dirichlet_bcs_quad& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord = default_ord);

double
nedelec_gradient (tmesh::quadrant_iterator & q,
                  const q1_vec& u, size_t i);

gradient
bim2c_quadtree_pde_recovered_gradient (tmesh& mesh,
                                       const q1_vec& u,
                                       active_fun is_active =
				       [] (tmesh::quadrant_iterator) {return true;});

q2_vec
bim2c_quadtree_pde_recovered_solution (tmesh& mesh,
                                       const q1_vec& u,
                                       const gradient& du);

double
estimator_grad (tmesh::quadrant_iterator q,
                const gradient & du_star,
                const q1_vec & u);

int 
zz_marker_grad (tmesh::quadrant_iterator q,
                const gradient & du_star,
                const q1_vec & u,
                double limit);

double
estimator_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const q1_vec & u);

int 
zz_marker_sol (tmesh::quadrant_iterator q,
               const q2_vec& ustar,
               const q1_vec& u,
               double limit);

double
l2_error (tmesh::quadrant_iterator q,
          const func & u_ex,
          const q1_vec & u);

double
semih1_error (tmesh::quadrant_iterator q,
              const func & dudx_ex,
              const func & dudy_ex,
              const q1_vec & u);

double
l2_star_error (tmesh::quadrant_iterator q,
               const func & u_ex,
               const q2_vec & ustar);

double
semih1_star_error (tmesh::quadrant_iterator q,
                   const func & dudx_ex,
                   const func & dudy_ex,
                   const gradient & du_star);

#endif
