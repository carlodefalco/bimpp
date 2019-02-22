#ifndef HAVE_QUAD_OPERATORS_3D_H
#define HAVE_QUAD_OPERATORS_3D_H 1

#include "bim_distributed_vector.h"
#include "bim_sparse.h"
#include "operators.h"
#include "tmesh_3d.h"

#include <functional>
#include <tuple>
#include <vector>

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
                 sparse_matrix& A);

void
bim3a_advection_diffusion (tmesh_3d & mesh,
                           const std::vector<double>& alpha,
                           const std::vector<double>& psi,
                           sparse_matrix& A);

void
bim3a_advection_eafe_diffusion (tmesh_3d & mesh,
                                const std::vector<double>& alpha,
                                const std::vector<double>& psi,
                                sparse_matrix& A);

void
bim3a_reaction (tmesh_3d& mesh,
                const std::vector<double>& delta,
                const std::vector<double>& zeta,
                sparse_matrix& A);

void
bim3a_rhs (tmesh_3d& mesh,
           const std::vector<double>& f,
           const std::vector<double>& g,
           std::vector<double>& rhs);

void
bim3a_solution_with_ghosts (tmesh_3d& mesh,
                            distributed_vector& rhs,
                            const binary_operator &op = std::plus<double> (),
                            const ordering& ord = default_ord);

std::vector<double>
bim3a_boundary_mass (tmesh_3d & mesh,
		     const int & tree_idx,
		     const int & boundary_idx,
		     std::vector<double> & M,
		     const func3_quad & fun =
		     [] (tmesh_3d::quadrant_iterator, tmesh_3d::idx_t)
		       {return 1;}
		     );

void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3& bcs,
                    sparse_matrix& A, std::vector<double>& rhs);

void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3_quad& bcs,
                    sparse_matrix& A, std::vector<double>& rhs);

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

template <class T>
double
l2_error (tmesh_3d::quadrant_iterator q,
          const func3 & u_ex,
          const T & u);

#endif
