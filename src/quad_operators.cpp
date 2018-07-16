#include "quad_operators.h"
#include "bim_distributed_vector.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <set>
#include <limits>
#include <iomanip>
#include <tuple>

// Compute harmonic mean of a and b.
static inline double
hm (const double& a, const double& b)
{ return 2 / (1 / a + 1 / b); }

// MPI_User_function.
static void replace(double *invec, double *inoutvec,
                    int *len, MPI_Datatype *dtype)
{
  for (int i = 0; i < *len; ++i)
    if (invec[i] != 0 && inoutvec[i] == 0)
      inoutvec[i] = invec[i];
}

static
std::array<std::array<double, 4>, 4> Aloc;

static
std::array<double, 4> rhsloc;

///  The assumed numbering for quadrant nodes is
///  the following :
///
///   2------------------3
///   |                  |
///   |                  |
///   |                  |
///   |                  |
///   0------------------1
///
///  Diagonal interactions are excluded via reduced
///  quadrature so only local elements
///
///     (0,1) (0,2) (1,3) (2,3) 
///
///  and their structurally symmetric contributions
///  are filled.

 
static void
assemble (tmesh::quadrant_iterator& quadrant,
          const std::array<std::array<double, 4>, 4>& locmat,
          sparse_matrix& A,
          const ordering& ordr = default_ord,
          const ordering& ordc = default_ord)
{

  std::vector<unsigned int> rows, cols;
  rows.reserve (2);
  cols.reserve (2);
  int i, j, r, c;
        
  for (i = 0; i < 4; ++i)
    {
      rows.clear ();
      if (! quadrant->is_hanging (i))
        rows.push_back (ordr (quadrant->gt (i)));
      else
        {
          rows.push_back (ordr (quadrant->gparent (0, i)));
          rows.push_back (ordr (quadrant->gparent (1, i)));
        }


      for (j = 0; j < 4; ++j)
        {
          if (j == 3 - i) continue;
          cols.clear ();
          if (! quadrant->is_hanging (j))
            cols.push_back (ordc (quadrant->gt (j)));
          else
            {
              cols.push_back (ordc (quadrant->gparent (0, j)));
              cols.push_back (ordc (quadrant->gparent (1, j)));
            }
              
          for (r = 0; r < rows.size (); ++r)
            for (c = 0; c < cols.size (); ++c)
              A[rows[r]][cols[c]] += locmat[i][j]
                / (rows.size () * cols.size ());
        }
    }

}


/// When matrix is lumped only diagonal entries need to
/// be assembled.
static void
assemble_diag (tmesh::quadrant_iterator& quadrant,
               const std::array<std::array<double, 4>, 4>& locmat,
               sparse_matrix& A,
               const ordering& ordr = default_ord,
               const ordering& ordc = default_ord)
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
        A[ordr (rows[r])][ordc (rows[r])] +=
          locmat[i][i] / rows.size ();
    }
}

/// Assemble rhs.
template <class T>
void
assemble_rhs (tmesh::quadrant_iterator& quadrant,
              const std::array<double, 4>& locrhs,
              T& rhs,
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

void
bim2a_structure (tmesh &tmsh,
                 sparse_matrix& A,
                 const ordering& ordr,
                 const ordering& ordc)
{ 
  for (auto row : Aloc)
    row.fill (0.0);
               
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    assemble (quadrant, Aloc, A, ordr, ordc);
  
  A.set_properties ();
}


void 
bim2a_advection_diffusion_loc
(tmesh::quadrant_iterator& quadrant,
 const double & alpha,
 const std::array<double, 4>& psi,
 std::array<std::array<double,4>,4>& locmat)
{
  double psi01 = 0;
  double psi13 = 0;
  double psi32 = 0;
  double psi20 = 0;
  
  double bp01 = 0, bm01 = 0;
  double bp13 = 0, bm13 = 0;
  double bp32 = 0, bm32 = 0;
  double bp20 = 0, bm20 = 0;
  
  double
    hx = quadrant->p (0, 1) - quadrant->p (0, 0),
    hy = quadrant->p (1, 2) - quadrant->p (1, 0);

  double
    hxby2hy = .5 * hx / hy,
    hyby2hx = .5 * hy / hx;
    
  std::vector<unsigned int> rows, cols;
  rows.reserve (2);
  cols.reserve (2);
  
  psi01 = psi[1] - psi[0];
  psi13 = psi[3] - psi[1];
  psi32 = psi[2] - psi[3];
  psi20 = psi[0] - psi[2];
  
  bimu_bernoulli (psi01, bp01, bm01);
  bimu_bernoulli (psi13, bp13, bm13);
  bimu_bernoulli (psi32, bp32, bm32);
  bimu_bernoulli (psi20, bp20, bm20);
  
  bp01 *= alpha * hyby2hx;
  bm01 *= alpha * hyby2hx;
  bp13 *= alpha * hxby2hy;
  bm13 *= alpha * hxby2hy;
  bp32 *= alpha * hyby2hx;
  bm32 *= alpha * hyby2hx;
  bp20 *= alpha * hxby2hy;
  bm20 *= alpha * hxby2hy;
  
  Aloc[0] = { bm01+bp20, -bp01,      -bm20,       0        };
  Aloc[1] = {-bm01,       bp01+bm13,  0,         -bp13     };
  Aloc[2] = {-bp20,       0,          bp32+bm20, -bm32     };
  Aloc[3] = { 0,         -bm13,      -bp32,       bm32+bp13};

}

template <class T>
void
bim2a_advection_diffusion (tmesh& mesh,
                           const std::vector<double>& alpha,
                           const T& psi,
                           sparse_matrix& A,                           
                           const ordering& ordr,
                           const ordering& ordc)
{
  for (auto row : Aloc)
    row.fill (0.0);
  
  double alpha_loc = 0;
  std::array<double, 4> psi_loc;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      alpha_loc = alpha[quadrant->get_forest_quad_idx ()];
      
      for (int n = 0; n < 4; ++n)
        {
          if (! quadrant->is_hanging (n))
            psi_loc[n] = psi[quadrant->gt (n)];
          else
            psi_loc[n] = 0.5 * (psi[quadrant->gparent (0, n)] +
                                psi[quadrant->gparent (1, n)]);
        }
  
      bim2a_advection_diffusion_loc (quadrant, alpha_loc, psi_loc, Aloc);
      assemble (quadrant, Aloc, A, ordr, ordc);
    }
}


void 
bim2a_advection_eafe_diffusion_loc
(tmesh::quadrant_iterator& quadrant,
 const std::array<double, 4>& alpha,
 const std::array<double, 4>& psi,
 std::array<std::array<double,4>,4>& locmat)
{
  double psi01 = 0;
  double psi13 = 0;
  double psi32 = 0;
  double psi20 = 0;
  
  double bp01 = 0, bm01 = 0;
  double bp13 = 0, bm13 = 0;
  double bp32 = 0, bm32 = 0;
  double bp20 = 0, bm20 = 0;
  
  double
    hx = quadrant->p (0, 1) - quadrant->p (0, 0),
    hy = quadrant->p (1, 2) - quadrant->p (1, 0);

  double
    hxby2hy = .5 * hx / hy,
    hyby2hx = .5 * hy / hx;
  
  psi01 = psi[1] - psi[0];
  psi13 = psi[3] - psi[1];
  psi32 = psi[2] - psi[3];
  psi20 = psi[0] - psi[2];
      
  bimu_bernoulli (psi01, bp01, bm01);
  bimu_bernoulli (psi13, bp13, bm13);
  bimu_bernoulli (psi32, bp32, bm32);
  bimu_bernoulli (psi20, bp20, bm20);
            
  bp01 *= hm (alpha[0], alpha[1]) * hyby2hx;
  bm01 *= hm (alpha[0], alpha[1]) * hyby2hx;
  bp13 *= hm (alpha[1], alpha[3]) * hxby2hy;
  bm13 *= hm (alpha[1], alpha[3]) * hxby2hy;
  bp32 *= hm (alpha[3], alpha[2]) * hyby2hx;
  bm32 *= hm (alpha[3], alpha[2]) * hyby2hx;
  bp20 *= hm (alpha[2], alpha[0]) * hxby2hy;
  bm20 *= hm (alpha[2], alpha[0]) * hxby2hy;
  
  locmat[0] = { bm01 + bp20, -bp01,        -bm20,         0          };
  locmat[1] = {-bm01,         bp01 + bm13,  0,           -bp13       };
  locmat[2] = {-bp20,         0,            bp32 + bm20, -bm32       };
  locmat[3] = { 0,           -bm13,        -bp32,         bm32 + bp13};
  
}


template <class T>
void 
bim2a_advection_eafe_diffusion (tmesh& mesh,
                                const T& alpha,
                                const T& psi,
                                sparse_matrix& A,
                                const ordering& ordr,
                                const ordering& ordc)
{
  for (auto row : Aloc)
    row.fill (0.0);
  
  std::array<double, 4> alpha_loc;
  std::array<double, 4> psi_loc;
      
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep (); ++quadrant)
    {
      for (int n = 0; n < 4; ++n)
        {
          if (! quadrant->is_hanging (n))
            {
              alpha_loc[n] = alpha[quadrant->gt (n)];
              psi_loc[n] = psi[quadrant->gt (n)];
            }
          else
            {
              alpha_loc[n] = 0.5 * (alpha[quadrant->gparent (0, n)] +
                                    alpha[quadrant->gparent (1, n)]);
              psi_loc[n] = 0.5 * (psi[quadrant->gparent (0, n)] +
                                  psi[quadrant->gparent (1, n)]);
            }
        }
      
      bim2a_advection_eafe_diffusion_loc (quadrant, alpha_loc, psi_loc, Aloc);
      assemble (quadrant, Aloc, A, ordr, ordc);      
    }
}

void
bim2a_reaction_loc (tmesh::quadrant_iterator& quadrant,
                    const double& delta,
                    const std::array<double, 4>& zeta,
                    std::array<std::array<double,4>,4>& locmat)
{
  auto hxhyby4 = .25 * (quadrant->p(0, 1) - quadrant->p (0, 0)) *
    (quadrant->p(1, 2) - quadrant->p (1, 0));
      
  for (int i = 0; i < 4; ++i)
    locmat[i][i] = (delta * zeta[i] * hxhyby4);
}

template <class T>
void
bim2a_reaction (tmesh& mesh,
                const std::vector<double>& delta,
                const T& zeta,
                sparse_matrix& A,
                const ordering& ordr,
                const ordering& ordc)
{  
  double delta_loc = 0;
  std::array<double, 4> zeta_loc;
  
  for (auto row : Aloc)
    row.fill (0.0);
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep (); ++quadrant)
    {
      delta_loc = delta[quadrant->get_forest_quad_idx ()];

      for (int n = 0; n < 4; ++n)
        {
          if (! quadrant->is_hanging (n))
            zeta_loc[n] = zeta[quadrant->gt (n)];
          else
            zeta_loc[n] = 0.5 * (zeta[quadrant->gparent (0, n)] +
                                 zeta[quadrant->gparent (1, n)]);
        }
      
      bim2a_reaction_loc (quadrant, delta_loc, zeta_loc, Aloc);        
      assemble_diag (quadrant, Aloc, A, ordr, ordc);
    }
}


void
bim2a_rhs_loc (tmesh::quadrant_iterator& quadrant,
               const double & f,
               const std::array<double, 4>& g,
               std::array<double, 4>& locrhs)
{
  auto hxhyby4 = .25 * (quadrant->p(0, 1) - quadrant->p (0, 0)) *
    (quadrant->p(1, 2) - quadrant->p (1, 0));
  
  for (int i = 0; i < 4; ++i)
    locrhs[i] = (f * g[i] * hxhyby4);
}

template <class T>
void
bim2a_rhs (tmesh& mesh,
           const std::vector<double>& f,
           const T& g,
           T& rhs,
           const ordering& ord)
{
  rhsloc.fill (0.0);
  
  double f_loc = 0;
  std::array<double, 4> g_loc;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep (); ++quadrant)
    {
      f_loc = f[quadrant->get_forest_quad_idx ()];

      for (int n = 0; n < 4; ++n)
        {
          if (! quadrant->is_hanging (n))
            g_loc[n] = g[quadrant->gt (n)];
          else
            g_loc[n] = 0.5 * (g[quadrant->gparent (0, n)] +
                              g[quadrant->gparent (1, n)]);
        }
      
      bim2a_rhs_loc (quadrant, f_loc, g_loc, rhsloc);
      assemble_rhs (quadrant, rhsloc, rhs, ord);
    }
}

template <class T>
void
bim2a_boundary_mass (tmesh& mesh,
                     const int & tree_idx,
                     const int & boundary_idx,
                     T & M,
                     const ordering & ord,
                     const func_quad & fun)
{
  double h = 0;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep (); ++quadrant)
    {
      if (quadrant->get_tree_idx () == tree_idx)
        {
          for (int i = 0; i < 4; ++i)
            {
              if (quadrant->e (i) == boundary_idx)
                {
                  if (boundary_idx == 0 ||
                      boundary_idx == 1)
                    h = quadrant->p(1, 2) - quadrant->p(1, 0);
                  else
                    h = quadrant->p(0, 1) - quadrant->p(0, 0);
                  
                  M[ord (quadrant->gt(i))] += 0.5 * h * fun (quadrant, i);
                }
            }
        }
    }
}


template <class T>
void
bim2a_dirichlet_bc_loc (sparse_matrix& A,
                        T& rhs,
                        const unsigned int& row,
                        const double& value,
                        const bool& only_rhs)
{
  rhs[row] = value;
  
  if (std::abs (A[row][row])
      < std::numeric_limits<double>::epsilon ())
    {
      A[row][row] = std::accumulate
        (A[row].begin (),
         A[row].end (),
         0.0,
         [] (double sum,
             const std::map<int, double>::value_type & p)
         {
           return (sum + std::abs (p.second));
         }
         );
    }

  if (! only_rhs)
    A[row][row] *= 1e16;
                    
  // Multiply rhs by the diagonal entry.
  rhs[row] *= A[row][row];
}


template <class T>
void
bim2a_dirichlet_bc (tmesh& mesh, const dirichlet_bcs& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord,
                    const bool& only_rhs)
{
  int boundary_idx, tree_idx;
  unsigned int row, col;
  
  std::set<unsigned int> marked;
  
  double value;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep (); ++quadrant)
    {
      tree_idx = quadrant->get_tree_idx ();
      
      for (int i = 0; i < 4; ++i)
        {
          boundary_idx = quadrant->e (i);
          row = ord (quadrant->gt (i));
          
          // If current node is on boundary and has not
          // been handled before.
          if (boundary_idx != tmesh::quadrant_t::NOT_ON_BOUNDARY
              && marked.count(row) == 0)
            {
              // Loop over all the boundary conditions.
              for (size_t bc = 0; bc < bcs.size (); ++bc)
                // If this boundary condition matches with
                // the current node.
                if (std::get<0> (bcs[bc]) == tree_idx
                    && std::get<1> (bcs[bc]) == boundary_idx)
                  {
                    // Mark current node so to avoid duplicate operations.
                    marked.insert (row);
                    
                    // Evaluate bc at current node.
                    value = (std::get<2> (bcs[bc]))
                      (quadrant->p (0, i), quadrant->p (1, i));
                    
                    bim2a_dirichlet_bc_loc (A, rhs, row, value, only_rhs);
                  }
            }
        }
    }
}

template <class T>
void
bim2a_dirichlet_bc (tmesh& mesh, const dirichlet_bcs_quad& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord,
                    const bool& only_rhs)
{
  int boundary_idx, tree_idx;
  unsigned int row, col;
  
  std::set<unsigned int> marked;

  double value;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep (); ++quadrant)
    {
      tree_idx = quadrant->get_tree_idx ();
      
      for (int i = 0; i < 4; ++i)
        {
          boundary_idx = quadrant->e (i);
          row = ord (quadrant->gt (i));
          
          // If current node is on boundary and has not
          // been handled before.
          if (boundary_idx != tmesh::quadrant_t::NOT_ON_BOUNDARY
              && marked.count(row) == 0)
            {
              // Loop over all the boundary conditions.
              for (size_t bc = 0; bc < bcs.size (); ++bc)
                // If this boundary condition matches with
                // the current node.
                if (std::get<0> (bcs[bc]) == tree_idx
                    && std::get<1> (bcs[bc]) == boundary_idx)
                  {
                    // Mark current node so to avoid duplicate operations.
                    marked.insert (row);
                    
                    // Evaluate bc at current node.
                    value = (std::get<2> (bcs[bc]))
                      (quadrant, i);
                    
                    bim2a_dirichlet_bc_loc (A, rhs, row, value, only_rhs);
                  }
            }
        }
    }
}

template <class T>
void
bim2a_robin_bc_loc (tmesh& mesh,
                    sparse_matrix& A,
                    T& rhs,
                    T& M_boundary,
                    const unsigned int& row,
                    const double& value_A,
                    const double& value_rhs)
{
  size_t start = mesh.lnodes->global_offset;
  size_t end = start + mesh.num_owned_nodes ();
  
  // If current node is owned.
  if (row >= start && row < end)
    rhs[row] = M_boundary[row] * value_rhs;
  else
    rhs[row] = 0;
              
  if (A[row].size ())
    {
      for (auto col = A[row].begin ();
           col != A[row].end ();
           ++col)
        {
          A[row][A.col_idx (col)] *= value_A;
        }
                  
      // If current node is owned.
      if (row >= start && row < end)
        A[row][row] += M_boundary[row];
    }
}

// Specialization.
template <>
void
interpolate_vector (tmesh & mesh,
                    std::vector<double> & vec_in,
                    std::vector<double> & vec_out,
                    const ordering & ord)
{
  tmesh::data_t * data;
  
  size_t start = mesh.lnodes->global_offset;
  size_t end = start + mesh.num_owned_nodes ();
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      data = static_cast<tmesh::data_t *> (quadrant->the_quadrant->p.user_data);
      
      for (int node = 0; node < 4; ++node)
        {
          // If current node is owned.
          if (! quadrant->is_hanging (node) &&
              vec_out[ord (quadrant->gt (node))] == 0 &&
              quadrant->gt (node) >= start && quadrant->gt (node) < end)
            {
              // Multiply by interpolation matrix.
              for (int i = 0; i < 4; ++i)
                vec_out[ord (quadrant->gt (node))] +=
                  data->interp_coeff[node][i] *
                  vec_in[ord (data->interp_idx[i])];
            }
        }
    }
}

// Specialization.
template <>
void
interpolate_vector (tmesh & mesh,
                    distributed_vector & vec_in,
                    distributed_vector & vec_out,
                    const ordering & ord)
{
  tmesh::data_t * data;
  
  // Assemble indices related to interpolation matrices.
  vec_in.clear_non_local ();
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      data = static_cast<tmesh::data_t *> (quadrant->the_quadrant->p.user_data);
      
      for (int node = 0; node < 4; ++node)
        {
          if (! quadrant->is_hanging (node))
            {
              // Multiply by interpolation matrix.
              for (int i = 0; i < 4; ++i)
                vec_in[ord (data->interp_idx[i])] += 0;
            }
        }
    }
  
  vec_in.assemble (replace_op);
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      data = static_cast<tmesh::data_t *> (quadrant->the_quadrant->p.user_data);
      
      for (int node = 0; node < 4; ++node)
        {
          if (! quadrant->is_hanging (node))
            {
              if (vec_out[ord (quadrant->gt (node))] == 0)
                {
                  // Multiply by interpolation matrix.
                  for (int i = 0; i < 4; ++i)
                    vec_out[ord (quadrant->gt (node))] +=
                      data->interp_coeff[node][i] *
                      vec_in[ord (data->interp_idx[i])];
                }
            }
          // Assemble parents.
          else
            {
              vec_out[ord (quadrant->gparent (0, node))] += 0;
              vec_out[ord (quadrant->gparent (1, node))] += 0;
            }
        }
    }
}

/// Edge ordering derived from vertex ordering
///
///   2----------------->3
///   ^        3         ^
///   |                  |
///   |0                1|
///   |                  |
///   |        2         |
///   0----------------->1

static constexpr
std::array<std::array<int, 3>, 4> edge =
  {0,2,1, 1,3,1, 0,1,0, 2,3,0};

template <class T>
static double
nedelec_gradient (tmesh::quadrant_iterator & q,
                  const T & u, size_t i)
{
  std::array<double, 4> u_aux;
  
  for (int n = 0; n < 4; ++n)
    {
      if (! q->is_hanging (n))
        u_aux[n] = u[q->gt (n)];
      else
        u_aux[n] = 0.5 * (u[q->gparent (0, n)] +
                          u[q->gparent (1, n)]);
    }

  double h = q->p (edge[i][2], edge[i][1]) -
    q->p (edge[i][2], edge[i][0]);
  
  return ((u_aux[edge[i][1]] - u_aux[edge[i][0]]) / h);
}

template <class T>
std::tuple<double, double, bool, bool>
bim2c_recovered_gradient_loc (tmesh::quadrant_iterator quadrant,
                              int node,
                              const T & u,
                              active_fun is_active)
{
  double hx = quadrant->p (0, 1) - quadrant->p (0, 0);
  double hy = quadrant->p (1, 2) - quadrant->p (1, 0);
  
  int node_n = 0, node_side = 0;
  std::vector<double> du_x, weights_x;
  std::vector<double> du_y, weights_y;
  
  du_x.clear (); weights_x.clear ();
  du_y.clear (); weights_y.clear ();

  double du_x_star = 0;
  double du_y_star = 0;
  
  bool assigned_x = false;
  bool assigned_y = false;
  
  // Compute Nedelec gradient on current element.
  switch (node)
    {
    case 0:
      du_x.push_back (nedelec_gradient (quadrant, u, 2));
      du_y.push_back (nedelec_gradient (quadrant, u, 0));
      break;
    case 1:
      du_x.push_back (nedelec_gradient (quadrant, u, 2));
      du_y.push_back (nedelec_gradient (quadrant, u, 1));
      break;
    case 2:
      du_x.push_back (nedelec_gradient (quadrant, u, 3));
      du_y.push_back (nedelec_gradient (quadrant, u, 0));
      break;
    case 3:
      du_x.push_back (nedelec_gradient (quadrant, u, 3));
      du_y.push_back (nedelec_gradient (quadrant, u, 1));
      break;
    }
          
  weights_x.push_back (1 / hx);
  weights_y.push_back (1 / hy);
  
  // Loop over face neighbors of current quadrant.
  for (auto neighbor = quadrant->begin_neighbor_sweep ();
       neighbor != quadrant->end_neighbor_sweep ();
       ++neighbor)
    {
      // Skip missing neighbors.
      if (neighbor->get_global_quad_idx () ==
          quadrant->get_global_quad_idx ())
        continue;
              
      // Skip inactive neighbors.
      if (! is_active(neighbor))
        continue;
              
      // Check if neighbor contains current vertex ("node").
      for (node_n = 0; node_n < 4; ++node_n)
        if (neighbor->gt (node_n) == quadrant->gt (node))
          break;
              
      // If not, or if node_n is hanging,
      // switch to the next neighbor.
      if (node_n == 4 || neighbor->is_hanging (node_n))
        continue;
              
      hx = neighbor->p (0, 1) - neighbor->p (0, 0);
      hy = neighbor->p (1, 2) - neighbor->p (1, 0);
              
      switch (node_n)
        {
        case 0:
          if (node == 1)
            du_x.push_back (nedelec_gradient(neighbor, u, 2));
          if (node == 2)
            du_y.push_back (nedelec_gradient(neighbor, u, 0));
          break;
        case 1:
          if (node == 0)
            du_x.push_back (nedelec_gradient(neighbor, u, 2));
          if (node == 3)
            du_y.push_back (nedelec_gradient(neighbor, u, 1));
          break;
        case 2:
          if (node == 3)
            du_x.push_back (nedelec_gradient(neighbor, u, 3));
          if (node == 0)
            du_y.push_back (nedelec_gradient(neighbor, u, 0));
          break;
        case 3:
          if (node == 2)
            du_x.push_back (nedelec_gradient(neighbor, u, 3));
          if (node == 1)
            du_y.push_back (nedelec_gradient(neighbor, u, 1));
          break;
        }
              
      if (weights_x.size () < du_x.size ())
        weights_x.push_back (1 / hx);
              
      if (weights_y.size () < du_y.size ())
        weights_y.push_back (1 / hy);
    }
  
  // If on any vertical boundary/interface.
  if (du_x.size () < 2)
    {
      for (auto neighbor = quadrant->begin_neighbor_sweep  ();
           neighbor != quadrant->end_neighbor_sweep ();
           ++neighbor)
        {
          // Skip missing neighbors.
          if (neighbor->get_global_quad_idx () ==
              quadrant->get_global_quad_idx ())
            continue;
                  
          // Skip inactive neighbors.
          if (! is_active (neighbor))
            continue;
                  
          // Check if neighbor contains the opposite vertex
          // of the horizontal side containing "node".
          switch (node)
            {
            case 0:
              node_side = 1;
              break;
            case 1:
              node_side = 0;
              break;
            case 2:
              node_side = 3;
              break;
            case 3:
              node_side = 2;
              break;
            }
                  
          for (node_n = 0; node_n < 4; ++node_n)
            if (neighbor->gt (node_n) ==
                quadrant->gt (node_side))
              break;
                  
          // If not, or if node_n is hanging,
          // switch to the next neighbor.
          if (node_n == 4 || neighbor->is_hanging (node_n))
            continue;
                  
          hx = neighbor->p (0, 1) - neighbor->p (0, 0);
                  
          switch (node_n)
            {
            case 0:
              if (node_side == 1)
                du_x.push_back (nedelec_gradient (neighbor, u, 2));
              break;
            case 1:
              if (node_side == 0)
                du_x.push_back (nedelec_gradient (neighbor, u, 2));
              break;
            case 2:
              if (node_side == 3)
                du_x.push_back (nedelec_gradient (neighbor, u, 3));
              break;
            case 3:
              if (node_side == 2)
                du_x.push_back (nedelec_gradient (neighbor, u, 3));
              break;
            }
                  
          if (weights_x.size () < du_x.size ())
            {
              weights_x[0] += 2 / hx;
              weights_x.push_back (-1 / hx);
            }
        }
    }
          
  // If on any horizontal boundary/interface.
  if (du_y.size () < 2)
    {
      for (auto neighbor = quadrant->begin_neighbor_sweep ();
           neighbor != quadrant->end_neighbor_sweep ();
           ++neighbor)
        {
          // Skip missing neighbors.
          if (neighbor->get_global_quad_idx () ==
              quadrant->get_global_quad_idx ())
            continue;
                  
          // Skip inactive neighbors.
          if (! is_active (neighbor))
            continue;
                  
          // Check if neighbor contains the opposite vertex
          // of the vertical side containing "node".
          switch (node)
            {
            case 0:
              node_side = 2;
              break;
            case 1:
              node_side = 3;
              break;
            case 2:
              node_side = 0;
              break;
            case 3:
              node_side = 1;
              break;
            }
                  
          for (node_n = 0; node_n < 4; ++node_n)
            if (neighbor->gt (node_n) ==
                quadrant->gt (node_side))
              break;
                  
          // If not, or if node_n is hanging,
          // switch to the next neighbor.
          if (node_n == 4 || neighbor->is_hanging (node_n))
            continue;
                  
          hy = neighbor->p (1, 2) - neighbor->p (1, 0);
                  
          switch (node_n)
            {
            case 0:
              if (node_side == 2)
                du_y.push_back (nedelec_gradient (neighbor, u, 0));
              break;
            case 1:
              if (node_side == 3)
                du_y.push_back (nedelec_gradient (neighbor, u, 1));
              break;
            case 2:
              if (node_side == 0)
                du_y.push_back (nedelec_gradient (neighbor, u, 0));
              break;
            case 3:
              if (node_side == 1)
                du_y.push_back (nedelec_gradient (neighbor, u, 1));
              break;
            }
                  
          if (weights_y.size () < du_y.size ())
            {
              weights_y[0] += 2 / hy;
              weights_y.push_back (-1 / hy);
            }
        }
    }
  
  assert (du_x.size () <= 2 && du_y.size () <= 2);
          
  if (du_x.size () == 2)
    {
      assigned_x = true;
      
      for (unsigned int ix = 0; ix < du_x.size (); ++ix)
        du_x_star +=
          du_x[ix] * weights_x[ix];
              
      du_x_star /=
        std::accumulate (weights_x.begin (),
                         weights_x.end (), 0.0);
    }
          
  if (du_y.size () == 2)
    {
      assigned_y = true;
      
      for (unsigned int iy = 0; iy < du_y.size (); ++iy)
        du_y_star +=
          du_y[iy] * weights_y[iy];
      
      du_y_star /=
        std::accumulate (weights_y.begin (),
                         weights_y.end (), 0.0);
    }
  
  return std::make_tuple (du_x_star, du_y_star,
                          assigned_x, assigned_y);
}

// Specialization.
template <>
gradient<std::vector<double>>
bim2c_quadtree_pde_recovered_gradient (tmesh & mesh,
                                       const std::vector<double> & u,
                                       active_fun is_active)
{
  std::vector<double> du_x_star (mesh.num_global_nodes (), 0);
  std::vector<double> du_y_star (mesh.num_global_nodes (), 0);  
  
  std::vector<bool> assigned_x (mesh.num_global_nodes (), 0);
  std::vector<bool> assigned_y (mesh.num_global_nodes (), 0);  
  
  std::tuple<double, double, bool, bool> du_star_loc;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      // Loop over non-hanging vertices of current quadrant.
      for (int node = 0; node < 4; ++node)
        {
          if (quadrant->is_hanging (node) ||
              (assigned_x[quadrant->gt (node)] &&
               assigned_y[quadrant->gt (node)]))
            continue;
          
          // Skip inactive quadrants.
          if (! is_active(quadrant))
            continue;
          
          du_star_loc =
            bim2c_recovered_gradient_loc (quadrant, node,
                                          u, is_active);
          
          du_x_star [quadrant->gt (node)] = std::get<0> (du_star_loc);
          du_y_star [quadrant->gt (node)] = std::get<1> (du_star_loc);
          
          assigned_x[quadrant->gt (node)] = std::get<2> (du_star_loc);
          assigned_y[quadrant->gt (node)] = std::get<3> (du_star_loc);
        }
    }
  
  // Send data to all processes so that non-assigned values
  // on current rank get assigned by other ranks.
  MPI_Op op;
  MPI_Op_create ((MPI_User_function *) replace, 1, &op);
  
  MPI_Allreduce (MPI_IN_PLACE, du_x_star.data (),
                 du_x_star.size (), MPI_DOUBLE,
                 op, MPI_COMM_WORLD);
  
  MPI_Allreduce (MPI_IN_PLACE, du_y_star.data (),
                 du_y_star.size (), MPI_DOUBLE,
                 op, MPI_COMM_WORLD);
  
  return std::make_pair (du_x_star, du_y_star);
}

// Specialization.
template <>
gradient<distributed_vector>
bim2c_quadtree_pde_recovered_gradient (tmesh & mesh,
                                       const distributed_vector & u,
                                       active_fun is_active)
{
  distributed_vector du_x_star (mesh.num_owned_nodes ());
  distributed_vector du_y_star (mesh.num_owned_nodes ());  
  
  distributed_vector assigned_x (mesh.num_owned_nodes ());
  distributed_vector assigned_y (mesh.num_owned_nodes ());  
  
  std::tuple<double, double, bool, bool> du_star_loc;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      // Loop over non-hanging vertices of current quadrant.
      for (int node = 0; node < 4; ++node)
        {
          if (assigned_x[quadrant->gt (node)] &&
              assigned_y[quadrant->gt (node)])
            continue;
          
          if (! quadrant->is_hanging (node))
            {
              // Assemble also entries related to inactive quadrants.
              du_x_star[quadrant->gt (node)] = 0;
              du_y_star[quadrant->gt (node)] = 0;
          
              // Skip inactive quadrants.
              if (! is_active(quadrant))
                continue;
          
              du_star_loc =
                bim2c_recovered_gradient_loc (quadrant, node,
                                              u, is_active);
          
              du_x_star [quadrant->gt (node)] = std::get<0> (du_star_loc);
              du_y_star [quadrant->gt (node)] = std::get<1> (du_star_loc);
          
              assigned_x[quadrant->gt (node)] = std::get<2> (du_star_loc);
              assigned_y[quadrant->gt (node)] = std::get<3> (du_star_loc);
            }
          // Assemble parents.
          else
            {
              du_x_star[quadrant->gparent (0, node)] += 0;
              du_x_star[quadrant->gparent (1, node)] += 0;
              
              du_y_star[quadrant->gparent (0, node)] += 0;
              du_y_star[quadrant->gparent (1, node)] += 0;
            }
        }
    }
  
  // Replace function for zero entries - i.e. those 
  // left unassigned on current process.
  binary_operator replace_zero =
    [] (const double & x, const double & y)
    {
      return (x != 0) ? x : y;
    };
  
  du_x_star.assemble (replace_zero);
  du_y_star.assemble (replace_zero);
  
  return std::make_pair (du_x_star, du_y_star);
}

template <class T>
q2_vec
bim2c_quadtree_pde_recovered_solution (tmesh & mesh,
                                       const T & u,
                                       const gradient<T> & du)
{
  q2_vec u_star (mesh.num_local_quadrants (),
                 std::array<double, 9> ({0,0,0,0,0,0,0,0,0}));
  
  double hx = 0, hy = 0;
  
  std::array<double, 4> u_star_loc,
    du_x_star_loc,
    du_y_star_loc;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p (0, 1) - quadrant->p (0, 0);
      hy = quadrant->p (1, 2) - quadrant->p (1, 0);
      
      // Compute values at vertices.
      for (int n = 0; n < 4; ++n)
        {
          if (! quadrant->is_hanging (n))
            {
              u_star_loc[n] = u[quadrant->gt (n)];
              du_x_star_loc[n] = du.first [quadrant->gt (n)];
              du_y_star_loc[n] = du.second[quadrant->gt (n)];
            }
          else
            {
              u_star_loc[n] = 0.5 * (u[quadrant->gparent (0, n)] +
                                     u[quadrant->gparent (1, n)]);
              
              // Determine whether n is hanging on an edge
              // directed along the x or y direction.
              int i = 0; int p = 0;
              
              for (; i < 4; ++i)
                {
                  if (quadrant->parent (0, n) == quadrant->t (i))
                    {
                      p = 0;
                      break;
                    }
                  else if (quadrant->parent (1, n) == quadrant->t (i))
                    {
                      p = 1;
                      break;
                    }
                }
              
              // Compute recovered solution at the
              // double-sized neighbor element.
              if (n == 0)
                {
                  if (i == 1)
                    u_star_loc[n] +=
                      (2 * hx) * (du.first[quadrant->gparent (1-p, n)] -
                                  du.first[quadrant->gparent (p, n)]) / 8;
                  else if (i == 2)
                    u_star_loc[n] +=
                      (2 * hy) * (du.second[quadrant->gparent (1-p, n)] -
                                  du.second[quadrant->gparent (p, n)]) / 8;
                }
              else if (n == 1)
                {
                  if (i == 0)
                    u_star_loc[n] +=
                      (2 * hx) * (du.first[quadrant->gparent (p, n)] -
                                  du.first[quadrant->gparent (1-p, n)]) / 8;
                  else if (i == 3)
                    u_star_loc[n] +=
                      (2 * hy) * (du.second[quadrant->gparent (1-p, n)] -
                                  du.second[quadrant->gparent (p, n)]) / 8;
                }
              else if (n == 2)
                {
                  if (i == 3)
                    u_star_loc[n] +=
                      (2 * hx) * (du.first[quadrant->gparent (1-p, n)] -
                                  du.first[quadrant->gparent (p, n)]) / 8;
                  else if (i == 0)
                    u_star_loc[n] +=
                      (2 * hy) * (du.second[quadrant->gparent (p, n)] -
                                  du.second[quadrant->gparent (1-p, n)]) / 8;
                }
              else if (n == 3)
                {
                  if (i == 2)
                    u_star_loc[n] +=
                      (2 * hx) * (du.first[quadrant->gparent (p, n)] -
                                  du.first[quadrant->gparent (1-p, n)]) / 8;
                  else if (i == 1)
                    u_star_loc[n] +=
                      (2 * hy) * (du.second[quadrant->gparent (p, n)] -
                                  du.second[quadrant->gparent (1-p, n)]) / 8;
                }
                
              du_x_star_loc[n] =
                0.5 * (du.first[quadrant->gparent (0, n)] +
                       du.first[quadrant->gparent (1, n)]);
              
              du_y_star_loc[n] =
                0.5 * (du.second[quadrant->gparent (0, n)] +
                       du.second[quadrant->gparent (1, n)]);
            }
          
          u_star[quadrant->get_forest_quad_idx ()][n] =
            u_star_loc[n];
        }
      
      // Compute values at faces.
      u_star[quadrant->get_forest_quad_idx ()][4] =
        0.5 * (u_star_loc[0] + u_star_loc[2])
        + hy * (du_y_star_loc[0] - du_y_star_loc[2]) / 8;
      
      u_star[quadrant->get_forest_quad_idx ()][5] =
        0.5 * (u_star_loc[1] + u_star_loc[3])
        + hy * (du_y_star_loc[1] - du_y_star_loc[3]) / 8;
      
      u_star[quadrant->get_forest_quad_idx ()][6] =
        0.5 * (u_star_loc[0] + u_star_loc[1])
        + hx * (du_x_star_loc[0] - du_x_star_loc[1]) / 8;
      
      u_star[quadrant->get_forest_quad_idx ()][7] =
        0.5 * (u_star_loc[2] + u_star_loc[3])
        + hx * (du_x_star_loc[2] - du_x_star_loc[3]) / 8;
      
      // Compute value at cell midpoint.
      u_star[quadrant->get_forest_quad_idx ()][8] =
        0.25 * (u_star[quadrant->get_forest_quad_idx ()][4] +
                u_star[quadrant->get_forest_quad_idx ()][5] +
                u_star[quadrant->get_forest_quad_idx ()][6] +
                u_star[quadrant->get_forest_quad_idx ()][7])
        + hx * 0.5 * (du_x_star_loc[0] + du_x_star_loc[2] -
                      du_x_star_loc[1] - du_x_star_loc[3]) / 16
        + hy * 0.5 * (du_y_star_loc[0] + du_y_star_loc[1] -
                      du_y_star_loc[2] - du_y_star_loc[3]) / 16;
    }
  
  return u_star;
}

// 4-points Gauss quadature nodes and weights (in [0, 1]).
static constexpr double gn[4] =
  {6.94318442029737e-02, 3.30009478207572e-01,
   6.69990521792428e-01, 9.30568155797026e-01};

static constexpr double gw[4] =
  {1.73927422568727e-01, 3.26072577431273e-01,
   3.26072577431273e-01, 1.73927422568727e-01};

// Transform nodes from [0, 1] to [x[0], x[1]].
static inline double
xformx (const double *x, const double X)
{ return (X * (x[1] - x[0]) + x[0]); }

// Transform weights from [0, 1] to [x[0], x[1]].
static inline double
xformw (const double *x, const double w)
{ return (w * (x[1] - x[0])); }

// Approximate integral of fun on [x[0], x[1]] x [y[0], y[1]].
static double
quad_integral (const double *x, const double *y,
               std::function<double (double, double)> fun)
{
  int ix, jy;
  double wx = 0, nx = 0, sum = 0;
  for (ix = 0; ix < 4; ++ix)
    {
      wx = xformw (x, gw[ix]);
      nx = xformx (x, gn[ix]);
      for (jy = 0; jy < 4; ++jy)        
        sum += fun (nx, xformx (y, gn[jy])) *
          wx * xformw (y, gw[jy]);
    }
  return (sum);
}

// Evaluate Nedelec x-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]])
// at (X, Y).
static double
dudx (double X, double Y, const double *x,
      const double *y, const double *u)
{
  double hxhy = (x[1] - x[0]) * (y[1] - y[0]);
  
  double db = (u[1] - u[0]);
  double dt = (u[3] - u[2]);
  
  return (db * (y[1] - Y) + dt * (Y - y[0])) / hxhy;
}

// Evaluate Nedelec y-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]])
// at (X, Y).
static double
dudy (double X, double Y, const double *x,
      const double *y, const double *u)
{
  double hxhy = (x[1] - x[0]) * (y[1] - y[0]);
  
  double dl = (u[2] - u[0]);
  double dr = (u[3] - u[1]);
  
  return (dl * (x[1] - X) + dr * (X - x[0])) / hxhy;
}

// Evaluate u (using Q1 basis functions
// on quadrant [x[0], x[1]] x [y[0], y[1]])
// at (X, Y).
static double
q1 (double X, double Y, const double *x,
    const double *y, const double *u)
{
  double hxhy = (x[1] - x[0]) * (y[1] - y[0]);
  double Xx1 = (X - x[1]);
  double Xx0 = (X - x[0]);
  double Yy1 = (Y - y[1]);
  double Yy0 = (Y - y[0]);
  return ((u[0] * Xx1  * Yy1 + u[1] * -Xx0 * Yy1 +
           u[2] * -Xx1 * Yy0 + u[3] * Xx0  * Yy0) / (hxhy));
}

// Evaluate u (using Q2 basis functions
// on quadrant [x[0], x[1]] x [y[0], y[1]])
// at (X, Y).
static double
q2 (double X, double Y, const double *x,
    const double *y, const double *u)
{
  double xc = 0.5 * (x[0] + x[1]);
  double yc = 0.5 * (y[0] + y[1]);

  double Xx1 = (X - x[1]);
  double Xx0 = (X - x[0]);
  double Xxc = (X - xc);
  double Yy1 = (Y - y[1]);
  double Yy0 = (Y - y[0]);
  double Yyc = (Y - yc);

  double hxhy2 = std::pow ((x[1] - x[0]) * (y[1] - y[0]), 2);

  return (u[0] * 4  * Xxc * Xx1 * Yyc * Yy1 +
          u[1] * 4  * Xx0 * Xxc * Yyc * Yy1 +
          u[2] * 4  * Xxc * Xx1 * Yy0 * Yyc +
          u[3] * 4  * Xx0 * Xxc * Yy0 * Yyc +
          u[4] * -8 * Xxc * Xx1 * Yy0 * Yy1 +
          u[5] * -8 * Xx0 * Xxc * Yy0 * Yy1 +
          u[6] * -8 * Xx0 * Xx1 * Yyc * Yy1 +
          u[7] * -8 * Xx0 * Xx1 * Yy0 * Yyc +
          u[8] * 16 * Xx0 * Xx1 * Yy0 * Yy1) / hxhy2;
}

// Compute ||grad^* u - grad u||_L^2(q).
template <class T>
double
estimator_grad (tmesh::quadrant_iterator q,
                const gradient<T> & du_star,
                const T & u)
{
  double
    x[2] = {q->p(0,0), q->p(0,1)},
    y[2] = {q->p(1,0), q->p(1,3)};
  
  double dudxstar_loc[4] = {0,0,0,0};
  double dudystar_loc[4] = {0,0,0,0};
  double u_loc[4] = {0,0,0,0};
  
  for (int ii = 0; ii < 4; ++ii)
    if (! q->is_hanging (ii))
      {
        dudxstar_loc[ii] = (du_star.first)[q->gt (ii)];
        dudystar_loc[ii] = (du_star.second)[q->gt (ii)];
        u_loc[ii] = u[q->gt(ii)];
      }
    else
      {
        dudxstar_loc[ii] = 0.5 *
          ((du_star.first)[q->gparent (0, ii)] +
           (du_star.first)[q->gparent (1, ii)]);
        
        dudystar_loc[ii] = 0.5 *
          ((du_star.second)[q->gparent (0, ii)] +
           (du_star.second)[q->gparent (1, ii)]);
        
        u_loc[ii] = 0.5 *
          (u[q->gparent (0, ii)] +
           u[q->gparent (1, ii)]);
      }

  
  auto fun =
    [x, y, dudxstar_loc, dudystar_loc, u_loc]
    (double X, double Y) -> double
    {
      return
      std::pow (dudx (X, Y, x, y, u_loc) -
                q1 (X, Y, x, y, dudxstar_loc), 2) +
      std::pow (dudy (X, Y, x, y, u_loc) -
                q1 (X, Y, x, y, dudystar_loc), 2);
    };
  
  return std::sqrt (quad_integral (x, y, fun));
}

// Refinement marker function based on ZZ estimator
// for the recovered gradient du*.
template <class T>
int
zz_marker_grad (tmesh::quadrant_iterator q,
                const gradient<T> & du_star,
                const T & u, double limit)
{ return estimator_grad (q, du_star, u) > limit ? 1 : 0; }

// Compute ||u^* - u||_L^2(q).
template <class T>
double
estimator_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const T & u)
{
  double
    x[2] = {q->p(0,0), q->p (0,1)},
    y[2] = {q->p(1,0), q->p (1,3)};

  double ustar_loc[9] = {0,0,0,0,0,0,0,0,0};
  double u_loc[4] = {0,0,0,0};

  for (int ii = 0; ii < 9; ++ii)
    ustar_loc[ii] = (ustar[q->get_forest_quad_idx ()])[ii];
  
  for (int ii = 0; ii < 4; ++ii)
    if (! q->is_hanging (ii))
      u_loc[ii] = u[q->gt (ii)];
    else
      u_loc[ii] = 0.5 *
        (u[q->gparent (0, ii)] +
         u[q->gparent (1, ii)]);

  auto fun =
    [x, y, ustar_loc, u_loc]
    (double X, double Y) -> double
    {
      return
      std::pow (q1 (X, Y, x, y, u_loc) -
                q2 (X, Y, x, y, ustar_loc), 2);
    };
    
  return std::sqrt (quad_integral (x, y, fun));
}

// Refinement marker function based on ZZ estimator
// for the recovered solution u*.
template <class T>
int
zz_marker_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const T & u,
               double limit)
{
  return estimator_sol (q, ustar, u) > limit ? 1 : 0;
}

// Compute ||u - u_ex||_L^2(q).
template <class T>
double
l2_error (tmesh::quadrant_iterator q,
          const func & u_ex,
          const T & u)
{
  double
    x[2] = {q->p (0,0), q->p (0,1)},
    y[2] = {q->p (1,0), q->p (1,3)};

  double u_loc[4] = {0,0,0,0};

  for (int ii = 0; ii < 4; ++ii)
    if (! q->is_hanging (ii))
      u_loc[ii] = u[q->gt (ii)];
    else
      u_loc[ii] = 0.5 *
        (u[q->gparent (0, ii)] + u[q->gparent (1, ii)]);

  auto fun =
    [x, y, u_loc, u_ex]
    (double X, double Y) -> double
    { return std::pow (q1 (X, Y, x, y, u_loc) - u_ex (X, Y), 2); };
    
  return std::sqrt (quad_integral (x, y, fun));
}

// Compute |u - u_ex|_H^1(q).
template <class T>
double
semih1_error (tmesh::quadrant_iterator q,
              const func & dudx_ex,
              const func & dudy_ex,
              const T & u)
{
  double
    x[2] = {q->p (0,0), q->p (0,1)},
    y[2] = {q->p (1,0), q->p (1,3)};

  double u_loc[4] = {0,0,0,0};

  for (int ii = 0; ii < 4; ++ii)
    if (! q->is_hanging (ii))
      u_loc[ii] = u[q->gt (ii)];
    else
      u_loc[ii] = 0.5 *
        (u[q->gparent (0, ii)] +
         u[q->gparent (1, ii)]);


  auto fun =
    [x, y, dudx_ex, dudy_ex, u_loc]
    (double X, double Y) -> double
    {
      return
      std::pow (dudx (X, Y, x, y, u_loc) -
                dudx_ex (X, Y), 2) +
      std::pow (dudy (X, Y, x, y, u_loc) -
                dudy_ex (X, Y), 2);
    };
    
  return std::sqrt (quad_integral (x, y, fun));
}

// Compute ||u_star - u_ex||_L^2(q).
double
l2_star_error (tmesh::quadrant_iterator q,
               const func & u_ex,
               const q2_vec & ustar)
{
  double
    x[2] = {q->p (0,0), q->p (0,1)},
    y[2] = {q->p (1,0), q->p (1,3)};
  
  double ustar_loc[9] = {0,0,0,0,0,0,0,0,0};
  
  for (int ii = 0; ii < 9; ++ii)
    ustar_loc[ii] = (ustar[q->get_forest_quad_idx ()])[ii];
  
  auto fun =
    [x, y, ustar_loc, u_ex]
    (double X, double Y) -> double
    {
      return
      std::pow (q2 (X, Y, x, y, ustar_loc) - u_ex (X, Y), 2);
    };
  
  return std::sqrt (quad_integral (x, y, fun));
}


// Compute ||du_star - grad(u_ex)||_L^2(q).
template <class T>
double
semih1_star_error (tmesh::quadrant_iterator q,
                   const func & dudx_ex,
                   const func & dudy_ex,
                   const gradient<T> & du_star)
{
  double
    x[2] = {q->p (0,0), q->p (0,1)},
    y[2] = {q->p (1,0), q->p (1,3)};
  
  double dudxstar_loc[4] = {0,0,0,0};
  double dudystar_loc[4] = {0,0,0,0};
  
  for (int ii = 0; ii < 4; ++ii)
    {
      if (! q->is_hanging (ii))
        {
          dudxstar_loc[ii] = (du_star.first)[q->gt (ii)];
          dudystar_loc[ii] = (du_star.second)[q->gt (ii)];
        }
      else
        {
          dudxstar_loc[ii] = 0.5 *
            ((du_star.first)[q->gparent (0, ii)] +
             (du_star.first)[q->gparent (1, ii)]);
          dudystar_loc[ii] = 0.5 *
            ((du_star.second)[q->gparent (0, ii)] +
             (du_star.second)[q->gparent (1, ii)]);
        }
    }
  
  auto fun =
    [x, y, dudxstar_loc, dudystar_loc, dudx_ex, dudy_ex]
    (double X, double Y) -> double
    {
      return
      std::pow (dudx_ex (X, Y) -
                q1 (X, Y, x, y, dudxstar_loc), 2) +
      std::pow (dudy_ex (X, Y) -
                q1 (X, Y, x, y, dudystar_loc), 2);
    };
  
  return std::sqrt (quad_integral (x, y, fun));
}

// Explicit instantiation.
template
void
assemble_rhs (tmesh::quadrant_iterator&,
              const std::array<double, 4>&,
              std::vector<double>&,
              const ordering& ord);

template
void
assemble_rhs (tmesh::quadrant_iterator&,
              const std::array<double, 4>&,
              distributed_vector&,
              const ordering& ord);

/* ---- */
template
void
bim2a_advection_diffusion (tmesh&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           sparse_matrix&,
                           const ordering&,
                           const ordering&);

template
void
bim2a_advection_diffusion (tmesh&,
                           const std::vector<double>&,
                           const distributed_vector&,
                           sparse_matrix&,
                           const ordering&,
                           const ordering&);

/* ---- */
template
void
bim2a_advection_eafe_diffusion (tmesh&,
                                const std::vector<double>&,
                                const std::vector<double>&,
                                sparse_matrix&,
                                const ordering&,
                                const ordering&);

template
void
bim2a_advection_eafe_diffusion (tmesh&,
                                const distributed_vector&,
                                const distributed_vector&,
                                sparse_matrix&,
                                const ordering&,
                                const ordering&);

/* ---- */
template
void
bim2a_reaction (tmesh&,
                const std::vector<double>&,
                const std::vector<double>&,
                sparse_matrix&,
                const ordering&,
                const ordering&);

template
void
bim2a_reaction (tmesh&,
                const std::vector<double>&,
                const distributed_vector&,
                sparse_matrix&,
                const ordering&,
                const ordering&);

/* ---- */
template
void
bim2a_rhs (tmesh&,
           const std::vector<double>&,
           const std::vector<double>&,
           std::vector<double>&,
           const ordering&);

template
void
bim2a_rhs (tmesh&,
           const std::vector<double>&,
           const distributed_vector&,
           distributed_vector&,
           const ordering&);

/* ---- */
template
void
bim2a_boundary_mass (tmesh&,
                     const int &,
                     const int &,
                     std::vector<double> &,
                     const ordering &,
                     const func_quad &);

template
void
bim2a_boundary_mass (tmesh&,
                     const int &,
                     const int &,
                     distributed_vector &,
                     const ordering &,
                     const func_quad &);

/* ---- */
template
void
bim2a_dirichlet_bc_loc (sparse_matrix&,
                        std::vector<double>&,
                        const unsigned int&,
                        const double&,
                        const bool&);

template
void
bim2a_dirichlet_bc_loc (sparse_matrix&,
                        distributed_vector&,
                        const unsigned int&,
                        const double&,
                        const bool&);

/* ---- */
template
void
bim2a_dirichlet_bc (tmesh&, const dirichlet_bcs&,
                    sparse_matrix&, std::vector<double>&,
                    const ordering&,
                    const bool&);

template
void
bim2a_dirichlet_bc (tmesh&, const dirichlet_bcs&,
                    sparse_matrix&, distributed_vector&,
                    const ordering&,
                    const bool&);

/* ---- */
template
void
bim2a_dirichlet_bc (tmesh&, const dirichlet_bcs_quad&,
                    sparse_matrix&, std::vector<double>&,
                    const ordering&,
                    const bool&);

template
void
bim2a_dirichlet_bc (tmesh&, const dirichlet_bcs_quad&,
                    sparse_matrix&, distributed_vector&,
                    const ordering&,
                    const bool&);

/* ---- */
template
void
bim2a_robin_bc_loc (tmesh& mesh,
                    sparse_matrix& A,
                    std::vector<double>& rhs,
                    std::vector<double>& M_boundary,
                    const unsigned int& row,
                    const double& value_A,
                    const double& value_rhs);

template
void
bim2a_robin_bc_loc (tmesh& mesh,
                    sparse_matrix& A,
                    distributed_vector& rhs,
                    distributed_vector& M_boundary,
                    const unsigned int& row,
                    const double& value_A,
                    const double& value_rhs);

/* ---- */
template
void
interpolate_vector (tmesh &,
                    std::vector<double> &,
                    std::vector<double> &,
                    const ordering &);

template
void
interpolate_vector (tmesh &,
                    distributed_vector &,
                    distributed_vector &,
                    const ordering &);

/* ---- */
template
double
nedelec_gradient (tmesh::quadrant_iterator &,
                  const std::vector<double> &, size_t);

template
double
nedelec_gradient (tmesh::quadrant_iterator &,
                  const distributed_vector &, size_t);

/* ---- */
template
std::tuple<double, double, bool, bool>
bim2c_recovered_gradient_loc (tmesh::quadrant_iterator,
                              int,
                              const std::vector<double> & u,
                              active_fun);

template
std::tuple<double, double, bool, bool>
bim2c_recovered_gradient_loc (tmesh::quadrant_iterator,
                              int,
                              const distributed_vector & u,
                              active_fun);

/* ---- */
template
gradient<std::vector<double>>
bim2c_quadtree_pde_recovered_gradient (tmesh &,
                                       const std::vector<double> &,
                                       active_fun);

template
gradient<distributed_vector>
bim2c_quadtree_pde_recovered_gradient (tmesh &,
                                       const distributed_vector &,
                                       active_fun);

/* ---- */
template
q2_vec
bim2c_quadtree_pde_recovered_solution (tmesh &,
                                       const std::vector<double> &,
                                       const gradient<std::vector<double>> &);
template
q2_vec
bim2c_quadtree_pde_recovered_solution (tmesh &,
                                       const distributed_vector &,
                                       const gradient<distributed_vector> &);
/* ---- */
template
double
estimator_grad (tmesh::quadrant_iterator q,
                const gradient<std::vector<double>> & du_star,
                const std::vector<double> & u);

template
double
estimator_grad (tmesh::quadrant_iterator q,
                const gradient<distributed_vector> & du_star,
                const distributed_vector & u);

/* ---- */
template
int
zz_marker_grad (tmesh::quadrant_iterator q,
                const gradient<std::vector<double>> & du_star,
                const std::vector<double> & u,
                double limit);

template
int
zz_marker_grad (tmesh::quadrant_iterator q,
                const gradient<distributed_vector> & du_star,
                const distributed_vector & u,
                double limit);

/* ---- */
template
double
estimator_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const std::vector<double> & u);

template
double
estimator_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const distributed_vector & u);

/* ---- */
template
int 
zz_marker_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const std::vector<double> & u,
               double limit);

template
int 
zz_marker_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const distributed_vector & u,
               double limit);

/* ---- */
template
double
l2_error (tmesh::quadrant_iterator q,
          const func & u_ex,
          const std::vector<double> & u);

template
double
l2_error (tmesh::quadrant_iterator q,
          const func & u_ex,
          const distributed_vector & u);

/* ---- */
template
double
semih1_error (tmesh::quadrant_iterator q,
              const func & dudx_ex,
              const func & dudy_ex,
              const std::vector<double> & u);

template
double
semih1_error (tmesh::quadrant_iterator q,
              const func & dudx_ex,
              const func & dudy_ex,
              const distributed_vector & u);

/* ---- */
template
double
semih1_star_error (tmesh::quadrant_iterator q,
                   const func & dudx_ex,
                   const func & dudy_ex,
                   const gradient<std::vector<double>> & du_star);
template
double
semih1_star_error (tmesh::quadrant_iterator q,
                   const func & dudx_ex,
                   const func & dudy_ex,
                   const gradient<distributed_vector> & du_star);
