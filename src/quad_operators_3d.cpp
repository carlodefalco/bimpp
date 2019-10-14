#include "quad_operators_3d.h"
#include <bim_timing.h>

#include <cmath>
#include <numeric>
#include <set>
#include <limits>
#include <iomanip>

static
std::array<std::array<double, 8>, 8> Aloc;

static void
assemble (tmesh_3d::quadrant_iterator& quadrant,
          const std::array<std::array<double, 8>, 8>& locmat,
          sparse_matrix& A,
          const ordering& ordr = default_ord,
          const ordering& ordc = default_ord)
{
  std::vector<unsigned int> rows, cols;
  std::vector<double> weights_r, weights_c;
  rows.reserve (4);
  cols.reserve (4);
  int i, j, r, c, pp;

  for  (i = 0; i < 8; ++i)
    {
      rows.clear ();

      if (! quadrant->is_hanging (i))
        rows.push_back (quadrant->gt (i));
      else
        for (pp = 0; pp < quadrant->num_parents (i); ++pp)
          rows.push_back (quadrant->gparent(pp, i));

      for(int j = 0; j < 8; ++j)
        {
          cols.clear ();

          if (! quadrant->is_hanging (j))
            cols.push_back (quadrant->gt (j));
          else
            for (pp = 0; pp < quadrant->num_parents (j); ++pp)
              cols.push_back (quadrant->gparent (pp, j));

          for (r = 0; r < rows.size (); ++r)
            if (locmat[i][j] != .0)
              for (int c = 0; c < cols.size (); ++c)
                {
                  A[rows[r]][cols[c]] += locmat[i][j] /
                    (rows.size () * cols.size ());
                }
        }
    }
}

// MPI_User_function.
static void replace (double *invec, double *inoutvec,
                     int *len, MPI_Datatype *dtype)
{
  for (int i = 0; i < *len; ++i)
    if (invec[i] != 0 && inoutvec[i] == 0)
      inoutvec[i] = invec[i];
}

void
bim3a_structure (tmesh_3d &tmsh,
                 sparse_matrix& A,
                 const ordering& ordr,
                 const ordering& ordc)
{
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    assemble (quadrant, Aloc, A, ordr, ordc);

  A.set_properties ();
}

template <class T>
void
bim3a_advection_diffusion (tmesh_3d& mesh,
                           const std::vector<double>& alpha,
                           const T& psi,
                           sparse_matrix& A,
                           bool symmetric,
                           const ordering& ordr,
                           const ordering& ordc)
{

  double bp01, bp13, bp23, bp02, bp04, bp15,
    bp26, bp37, bp45, bp57, bp67 ,bp46;

  double bm01, bm13, bm23, bm02, bm04, bm15,
    bm26, bm37, bm45, bm57, bm67 ,bm46;

  double hx, hy, hz;

  unsigned int iel = 0;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      std::array<double, 8> psi_aux;

      for (int n = 0; n < 8; ++n)
        {
          psi_aux[n] = 0;
          if (! quadrant->is_hanging (n))
            psi_aux[n] = psi[quadrant->gt (n)];
          else
            for (int pp = 0; pp < quadrant->num_parents (n); ++pp)
              psi_aux[n] += psi[quadrant->gparent (pp, n)] /
                quadrant->num_parents (n);
        }

      bimu_bernoulli(psi_aux[1] - psi_aux[0], bp01, bm01);
      bimu_bernoulli(psi_aux[3] - psi_aux[1], bp13, bm13);
      bimu_bernoulli(psi_aux[3] - psi_aux[2], bp23, bm23);
      bimu_bernoulli(psi_aux[2] - psi_aux[0], bp02, bm02);
      bimu_bernoulli(psi_aux[4] - psi_aux[0], bp04, bm04);
      bimu_bernoulli(psi_aux[5] - psi_aux[1], bp15, bm15);
      bimu_bernoulli(psi_aux[6] - psi_aux[2], bp26, bm26);
      bimu_bernoulli(psi_aux[7] - psi_aux[3], bp37, bm37);
      bimu_bernoulli(psi_aux[5] - psi_aux[4], bp45, bm45);
      bimu_bernoulli(psi_aux[7] - psi_aux[5], bp57, bm57);
      bimu_bernoulli(psi_aux[7] - psi_aux[6], bp67, bm67);
      bimu_bernoulli(psi_aux[6] - psi_aux[4], bp46, bm46);

      hx = quadrant->p(0, 7) - quadrant->p(0, 0);
      hy = quadrant->p(1, 7) - quadrant->p(1, 0);
      hz = quadrant->p(2, 7) - quadrant->p(2, 0);

      iel = quadrant->get_forest_quad_idx();
      double l01 = alpha[iel] * hy * hz / (4 * hx);
      double l02 = alpha[iel] * hx * hz / (4 * hy);
      double l04 = alpha[iel] * hx * hy / (4 * hz);
      double l13 = l02;
      double l15 = l04;
      double l23 = l01;
      double l26 = l04;
      double l37 = l04;
      double l45 = l01;
      double l46 = l02;
      double l57 = l02;
      double l67 = l01;

      bp01 *= l01; bm01 *= l01;
      bp02 *= l02; bm02 *= l02;
      bp04 *= l04; bm04 *= l04;
      bp13 *= l13; bm13 *= l13;
      bp15 *= l15; bm15 *= l15;
      bp23 *= l23; bm23 *= l23;
      bp26 *= l26; bm26 *= l26;
      bp37 *= l37; bm37 *= l37;
      bp45 *= l45; bm45 *= l45;
      bp46 *= l46; bm46 *= l46;
      bp57 *= l57; bm57 *= l57;
      bp67 *= l67; bm67 *= l67;

      Aloc[0] = {bm01+bm02+bm04, -bp01,      -bp02,      0.,
                 -bp04,      0.,         0.,         0.         };
      Aloc[1] = {    -bm01,  bp01+bm13+bm15, 0.,         -bp13,
                     0.,         -bp15,      0.,         0.         };
      Aloc[2] = {    -bm02,      0.,     bp02+bm23+bm26, -bp23,
                     0.,         0.,         -bp26,      0.         };
      Aloc[3] = {    0.,         -bm13,      -bm23,  bp13+bp23+bm37,
                     0.,         0.,         0.,         -bp37      };
      Aloc[4] = {    -bm04,      0.,         0.,         0.,
                     bp04+bm45+bm46, -bp45,      -bp46,      0.         };
      Aloc[5] = {    0.,         -bm15,      0.,         0.,
                     -bm45,  bp15+bp45+bm57, 0.,         -bp57      };
      Aloc[6] = {    0.,         0.,         -bm26,      0.,
                     -bm46,      0.,     bp26+bp46+bm67, -bp67      };
      Aloc[7] = {    0.,         0.,         0.,         -bm37,
                     0.,         -bm57,      -bm67,  bp37+bp57+bp67 };

      assemble (quadrant, Aloc, A, ordr, ordc);
    }
}

template <class T>
void
bim3a_reaction (tmesh_3d& mesh,
                const std::vector<double>& delta,
                const T& zeta,
                sparse_matrix& A,
                const ordering& ordr,
                const ordering& ordc)
{
  double hx, hy, hz;

  unsigned int iel = 0;
  std::vector<unsigned int> rows;
  rows.reserve(4);

  double z_loc = 0;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p(0, 7) - quadrant->p(0, 0);
      hy = quadrant->p(1, 7) - quadrant->p(1, 0);
      hz = quadrant->p(2, 7) - quadrant->p(2, 0);

      iel = quadrant->get_forest_quad_idx ();

      for(int i = 0; i < 8; ++i)
        {
          rows.clear();
          z_loc = 0;
          if (!quadrant->is_hanging (i))
            {
              rows.push_back (quadrant->gt (i));
              z_loc = zeta[quadrant->gt (i)];
            }
          else
            for (int pp = 0; pp < quadrant->num_parents (i); ++pp)
              {
                rows.push_back (quadrant->gparent (pp, i));
                z_loc += zeta[quadrant->gparent (pp, i)] /
                  quadrant->num_parents (i);
              }

          for (int r = 0; r < rows.size (); ++r)
            if (delta[iel] * z_loc != .0 )
              A[rows[r]][rows[r]] +=
                (delta[iel] * z_loc * hx * hy * hz / 8) / rows.size ();
        }
    }
}

template <class T>
void
bim3a_rhs (tmesh_3d& mesh,
           const std::vector<double>& f,
           const T& g, T& rhs,
           const ordering& ord)
{
  double hx, hy, hz;

  unsigned int iel = 0;
  std::vector<unsigned int> rows;
  rows.reserve (4);

  double g_loc = 0;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p (0, 7) - quadrant->p (0, 0);
      hy = quadrant->p (1, 7) - quadrant->p (1, 0);
      hz = quadrant->p (2, 7) - quadrant->p (2, 0);

      iel = quadrant->get_forest_quad_idx ();

      for(int i = 0; i < 8; ++i)
        {
          rows.clear ();
          g_loc = 0;
          if (! quadrant->is_hanging (i))
            {
              rows.push_back (quadrant->gt (i));
              g_loc = g[quadrant->gt (i)];
            }
          else
            for (int pp = 0; pp < quadrant->num_parents (i); ++pp)
              {
                rows.push_back (quadrant->gparent (pp, i));
                g_loc += g[quadrant->gparent (pp, i)] /
                  quadrant->num_parents (i);
              }

          for (int r = 0; r < rows.size(); ++r)
            rhs[rows[r]] +=
              (f[iel] * g_loc * hx * hy *hz / 8) / rows.size ();
        }
    }
}

void
bim3a_solution_with_ghosts (tmesh_3d& mesh,
                            distributed_vector& v,
                            const binary_operator &op,
                            const ordering& ord,
                            bool ra)
{
  int node = 0;
  for (auto q = mesh.begin_quadrant_sweep ();
       q != mesh.end_quadrant_sweep ();
       ++q)
    {
      for (node = 0; node < 8; ++node)
        {
          if (! q->is_hanging (node))
            v[ord(q->gt (node))] += 0;
        }

      for (auto n = q->begin_neighbor_sweep ();
           n != q->end_neighbor_sweep ();
           ++n)
        for (node = 0; node < 8; ++node)
          if (! n->is_hanging (node))
            v[ord(n->gt (node))] += 0;
          else
            {
              int np = n->num_parents(node);
              for (int pp = 0; pp < np; ++pp)
                v[ord(n->gparent(pp,node))] += 0;
            }
    }
  if (ra)
    {
      v.remap ();
      v.assemble (op);
    }
}

std::vector<double>
bim3a_boundary_mass (tmesh_3d & mesh,
                     const int & tree_idx,
                     const int & boundary_idx,
                     std::vector<double> & M,
                     const func3_quad & fun,
                     const ordering & ord)
{
  double area = 0;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      if (quadrant->get_tree_idx () == tree_idx)
        {
          for (int i = 0; i < 8; ++i)
            {
              if (quadrant->e (i) == boundary_idx)
                {
                  if (boundary_idx == 0 ||
                      boundary_idx == 1)
                    area =
                      (quadrant->p(1, 2) - quadrant->p(1, 0)) *
                      (quadrant->p(2, 4) - quadrant->p(2, 0));
                  else if (boundary_idx == 2 ||
                           boundary_idx == 3)
                    area =
                      (quadrant->p(0, 1) - quadrant->p(0, 0)) *
                      (quadrant->p(2, 4) - quadrant->p(2, 0));
                  else
                    area =
                      (quadrant->p(0, 1) - quadrant->p(0, 0)) *
                      (quadrant->p(1, 2) - quadrant->p(1, 0));

                  M[quadrant->gt(i)] += 0.25 * area * fun (quadrant, i);
                }
            }
        }
    }

  return M;
}

template <class T>
void
bim3a_dirichlet_bc_loc (sparse_matrix& A,
                        T& rhs,
                        const unsigned int& row,
                        const double& value,
                        const bool& only_rhs)
{
  if (std::abs (A[row][row]) < std::numeric_limits<double>::epsilon())
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

  // Multiply rhs by the diagonal entry
  rhs[row] = A[row][row] * value;
}

template <class T>
void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord,
                    const bool& only_rhs)
{
  int boundary_idx, tree_idx;
  unsigned int row, col;

  std::set<unsigned int> marked;

  double value;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      tree_idx = quadrant->get_tree_idx ();

      for (int i = 0; i < 8; ++i)
        {
          boundary_idx = quadrant->e (i);
          row = quadrant->gt (i);

          // If current node is on boundary and has not
          // been handled before.
          if (boundary_idx != tmesh_3d::quadrant_t::NOT_ON_BOUNDARY
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

                    // Evaluate bc at current node
                    value = (std::get<2> (bcs[bc]))
                            (quadrant->p (0, i),
                              quadrant->p (1, i),
                              quadrant->p (2, i));

                    bim3a_dirichlet_bc_loc (A, rhs, row, value, only_rhs);
                  }
            }
        }
    }
}

template <class T>
void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3_quad& bcs,
                    sparse_matrix& A, T& rhs,
                    const ordering& ord,
                    const bool& only_rhs)
{
  int boundary_idx, tree_idx;
  unsigned int row, col;

  std::set<unsigned int> marked;

  double value;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      tree_idx = quadrant->get_tree_idx ();

      for (int i = 0; i < 8; ++i)
        {
          boundary_idx = quadrant->e (i);
          row = quadrant->gt (i);

          // If current node is on boundary and has not
          // been handled before.
          if (boundary_idx != tmesh_3d::quadrant_t::NOT_ON_BOUNDARY
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

                    // Evaluate bc at current node
                    value =
                      (std::get<2> (bcs[bc])) (quadrant, i);

                    bim3a_dirichlet_bc_loc (A, rhs, row, value, only_rhs);                    
                  }
            }
        }
    }
}

// Specialization.
template <>
void
interpolate_vector (tmesh_3d & mesh,
                    std::vector<double> & vec_in,
                    std::vector<double> & vec_out,
                    const ordering & ord)
{
  tmesh_3d::data_t * data;

  size_t start = mesh.lnodes->global_offset;
  size_t end = start + mesh.num_owned_nodes ();

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      data = static_cast<tmesh_3d::data_t *> (quadrant->the_quadrant->p.user_data);

      for (int node = 0; node < 8; ++node)
        {
          // If current node is owned.
          if (! quadrant->is_hanging (node) &&
              vec_out[ord (quadrant->gt (node))] == 0 &&
              quadrant->gt (node) >= start && quadrant->gt (node) < end)
            {
              // Multiply by interpolation matrix.
              for (int i = 0; i < 8; ++i)
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
interpolate_vector (tmesh_3d & mesh,
                    distributed_vector & vec_in,
                    distributed_vector & vec_out,
                    const ordering & ord)
{
  tmesh_3d::data_t * data;

  // Assemble indices related to interpolation matrices.
  vec_in.clear_non_local ();

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      data = static_cast<tmesh_3d::data_t *> (quadrant->the_quadrant->p.user_data);

      for (int node = 0; node < 8; ++node)
        {
          if (! quadrant->is_hanging (node))
            {
              // Multiply by interpolation matrix.
              for (int i = 0; i < 8; ++i)
                vec_in[ord (data->interp_idx[i])] += 0;
            }
        }
    }

  vec_in.assemble (replace_op);

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      data = static_cast<tmesh_3d::data_t *> (quadrant->the_quadrant->p.user_data);

      for (int node = 0; node < 8; ++node)
        {
          if (! quadrant->is_hanging (node))
            {
              if (vec_out[ord (quadrant->gt (node))] == 0)
                {
                  // Multiply by interpolation matrix.
                  for (int i = 0; i < 8; ++i)
                    vec_out[ord (quadrant->gt (node))] +=
                      data->interp_coeff[node][i] *
                      vec_in[ord (data->interp_idx[i])];
                }
            }
          // Assemble parents.
          else
            {
              for (int pp = 0; pp < quadrant->num_parents(node); ++pp)
                vec_out[ord (quadrant->gparent (pp, node))] += 0;
            }
        }
    }
}

/// Edge ordering:
///
///             ______11_________
///            /|                /|
///           3 |               7 |
///          /  |              /  |
///         /___|____9________/   |
///        |    |             |   5
///        |    0             |   |
///        |    |             4   |
///        1    |_______10____|___|
///        |   /              |   /
///        |  2               |  6
///        | /                | /
///        |/_______8_________|/
///
///

static constexpr
std::array<std::array<int, 3>, 12> edge = {2,6,2, 0,4,2, 0,2,1, 4,6,1,
                                           1,5,2, 3,7,2, 1,3,1, 5,7,1,
                                           0,1,0, 4,5,0, 2,3,0, 6,7,0};

template <class T>
static double
nedelec_gradient (tmesh_3d::quadrant_iterator & q,
                  const T& u, size_t i)
{
  std::array<double,8> u_aux;

  for (int n = 0; n < 8; ++n)
    {
      if (! q->is_hanging(n))
        u_aux[n] = u[q -> gt(n)];
      else
        {
          u_aux[n] = 0.;
          int np = q->num_parents(n);
          for (int pp = 0; pp < np; ++pp)
            u_aux[n] += u[q->gparent(pp,n)];
          u_aux[n] /= np;
        }
    }

  double h = q->p (edge[i][2], edge[i][1]) -
    q->p (edge[i][2], edge[i][0]);

  return ((u_aux[edge[i][1]] - u_aux[edge[i][0]]) / h);
}

template <class T>
std::tuple<double, double, double, bool, bool, bool>
bim3c_recovered_gradient_loc (tmesh_3d::quadrant_iterator quadrant,
                              int node,
                              const T& u,
                              active_fun3 is_active)
{
  double hx = quadrant->p (0, 7) - quadrant->p (0, 0);
  double hy = quadrant->p (1, 7) - quadrant->p (1, 0);
  double hz = quadrant->p (2, 7) - quadrant->p (2, 0);

  int node_n = 0, node_side = 0;
  std::vector<double> du_x, weights_x;
  std::vector<double> du_y, weights_y;
  std::vector<double> du_z, weights_z;

  du_x.clear (); weights_x.clear ();
  du_y.clear (); weights_y.clear ();
  du_z.clear (); weights_z.clear ();

  du_x.reserve(2); weights_x.reserve(2);
  du_y.reserve(2); weights_y.reserve(2);
  du_z.reserve(2); weights_z.reserve(2);

  double du_x_star = 0;
  double du_y_star = 0;
  double du_z_star = 0;

  bool assigned_x = false;
  bool assigned_y = false;
  bool assigned_z = false;

  // Compute Nedelec gradient on current element.

  switch (node)
    {
    case 0:
      du_x.push_back (nedelec_gradient (quadrant, u, 8));
      du_y.push_back (nedelec_gradient (quadrant, u, 2));
      du_z.push_back (nedelec_gradient (quadrant, u, 1));
      break;
    case 1:
      du_x.push_back (nedelec_gradient (quadrant, u, 8));
      du_y.push_back (nedelec_gradient (quadrant, u, 6));
      du_z.push_back (nedelec_gradient (quadrant, u, 4));
      break;
    case 2:
      du_x.push_back (nedelec_gradient (quadrant, u, 10));
      du_y.push_back (nedelec_gradient (quadrant, u, 2));
      du_z.push_back (nedelec_gradient (quadrant, u, 0));
      break;
    case 3:
      du_x.push_back (nedelec_gradient (quadrant, u, 10));
      du_y.push_back (nedelec_gradient (quadrant, u, 6));
      du_z.push_back (nedelec_gradient (quadrant, u, 5));
      break;
    case 4:
      du_x.push_back (nedelec_gradient (quadrant, u, 9));
      du_y.push_back (nedelec_gradient (quadrant, u, 3));
      du_z.push_back (nedelec_gradient (quadrant, u, 1));
      break;
    case 5:
      du_x.push_back (nedelec_gradient (quadrant, u, 9));
      du_y.push_back (nedelec_gradient (quadrant, u, 7));
      du_z.push_back (nedelec_gradient (quadrant, u, 4));
      break;
    case 6:
      du_x.push_back (nedelec_gradient (quadrant, u, 11));
      du_y.push_back (nedelec_gradient (quadrant, u, 3));
      du_z.push_back (nedelec_gradient (quadrant, u, 0));
      break;
    case 7:
      du_x.push_back (nedelec_gradient (quadrant, u, 11));
      du_y.push_back (nedelec_gradient (quadrant, u, 7));
      du_z.push_back (nedelec_gradient (quadrant, u, 5));
      break;
    }

  weights_x.push_back (1 / hx);
  weights_y.push_back (1 / hy);
  weights_z.push_back (1 / hz);

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
      for (node_n = 0; node_n < 8; ++node_n)
        if (neighbor->gt (node_n) == quadrant->gt (node))
          break;

      // If not, or if node_n is hanging,
      // switch to the next neighbor.
      if (node_n == 8 || neighbor->is_hanging (node_n))
        continue;

      hx = neighbor->p (0, 7) - neighbor->p (0, 0);
      hy = neighbor->p (1, 7) - neighbor->p (1, 0);
      hz = neighbor->p (2, 7) - neighbor->p (2, 0);

      switch (node_n)
        {
        case 0:
          if (node == 1)
            du_x.push_back (nedelec_gradient(neighbor, u, 8));
          if (node == 2)
            du_y.push_back (nedelec_gradient(neighbor, u, 2));
          if (node == 4)
            du_z.push_back (nedelec_gradient(neighbor, u, 1));
          break;
        case 1:
          if (node == 0)
            du_x.push_back (nedelec_gradient(neighbor, u, 8));
          if (node == 3)
            du_y.push_back (nedelec_gradient(neighbor, u, 6));
          if (node == 5)
            du_z.push_back (nedelec_gradient(neighbor, u, 4));
          break;
        case 2:
          if (node == 3)
            du_x.push_back (nedelec_gradient(neighbor, u, 10));
          if (node == 0)
            du_y.push_back (nedelec_gradient(neighbor, u, 2));
          if (node == 6)
            du_z.push_back (nedelec_gradient(neighbor, u, 0));
          break;
        case 3:
          if (node == 2)
            du_x.push_back (nedelec_gradient(neighbor, u, 10));
          if (node == 1)
            du_y.push_back (nedelec_gradient(neighbor, u, 6));
          if (node == 7)
            du_z.push_back (nedelec_gradient(neighbor, u, 5));
          break;
        case 4:
          if (node == 5)
            du_x.push_back (nedelec_gradient(neighbor, u, 9));
          if (node == 6)
            du_y.push_back (nedelec_gradient(neighbor, u, 3));
          if (node == 0)
            du_z.push_back (nedelec_gradient(neighbor, u, 1));
          break;
        case 5:
          if (node == 4)
            du_x.push_back (nedelec_gradient(neighbor, u, 9));
          if (node == 7)
            du_y.push_back (nedelec_gradient(neighbor, u, 7));
          if (node == 1)
            du_z.push_back (nedelec_gradient(neighbor, u, 4));
          break;
        case 6:
          if (node == 7)
            du_x.push_back (nedelec_gradient(neighbor, u, 11));
          if (node == 4)
            du_y.push_back (nedelec_gradient(neighbor, u, 3));
          if (node == 2)
            du_z.push_back (nedelec_gradient(neighbor, u, 0));
          break;
        case 7:
          if (node == 6)
            du_x.push_back (nedelec_gradient(neighbor, u, 11));
          if (node == 5)
            du_y.push_back (nedelec_gradient(neighbor, u, 7));
          if (node == 3)
            du_z.push_back (nedelec_gradient(neighbor, u, 5));
          break;
        }

      if (weights_x.size () < du_x.size ())
        weights_x.push_back (1 / hx);

      if (weights_y.size () < du_y.size ())
        weights_y.push_back (1 / hy);

      if (weights_z.size () < du_z.size ())
        weights_z.push_back (1 / hz);
    }

  // If on any boundary/interface normal to x
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
          // of the side along x containing "node".
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
            case 4:
              node_side = 5;
              break;
            case 5:
              node_side = 4;
              break;
            case 6:
              node_side = 7;
              break;
            case 7:
              node_side = 6;
              break;
            }

          if (neighbor->gt (node) != quadrant->gt (node_side) ||
              neighbor->is_hanging (node))
            continue;

          hx = neighbor->p (0, 1) - neighbor->p (0, 0);

          switch (node)
            {
            case 0:
              if (node_side == 1)
                du_x.push_back (nedelec_gradient (neighbor, u, 8));
              break;
            case 1:
              if (node_side == 0)
                du_x.push_back (nedelec_gradient (neighbor, u, 8));
              break;
            case 2:
              if (node_side == 3)
                du_x.push_back (nedelec_gradient (neighbor, u, 10));
              break;
            case 3:
              if (node_side == 2)
                du_x.push_back (nedelec_gradient (neighbor, u, 10));
              break;
            case 4:
              if (node_side == 5)
                du_x.push_back (nedelec_gradient (neighbor, u, 9));
              break;
            case 5:
              if (node_side == 4)
                du_x.push_back (nedelec_gradient (neighbor, u, 9));
              break;
            case 6:
              if (node_side == 7)
                du_x.push_back (nedelec_gradient (neighbor, u, 11));
              break;
            case 7:
              if (node_side == 6)
                du_x.push_back (nedelec_gradient (neighbor, u, 11));
              break;
            }

          if (weights_x.size () < du_x.size ())
            {
              weights_x[0] += 2 / hx;
              weights_x.push_back (-1 / hx);
            }
        }
    }

  // If on any boundary/interface normal to y
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
          // of the side along y containing "node".
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
            case 4:
              node_side = 6;
              break;
            case 5:
              node_side = 7;
              break;
            case 6:
              node_side = 4;
              break;
            case 7:
              node_side = 5;
              break;
            }

          if (neighbor->gt (node) != quadrant->gt (node_side) ||
              neighbor->is_hanging (node))
            continue;

          hy = neighbor->p (1, 2) - neighbor->p (1, 0);

          switch (node)
            {
            case 0:
              if (node_side == 2)
                du_y.push_back (nedelec_gradient (neighbor, u, 2));
              break;
            case 1:
              if (node_side == 3)
                du_y.push_back (nedelec_gradient (neighbor, u, 6));
              break;
            case 2:
              if (node_side == 0)
                du_y.push_back (nedelec_gradient (neighbor, u, 2));
              break;
            case 3:
              if (node_side == 1)
                du_y.push_back (nedelec_gradient (neighbor, u, 6));
              break;
            case 4:
              if (node_side == 6)
                du_y.push_back (nedelec_gradient (neighbor, u, 3));
              break;
            case 5:
              if (node_side == 7)
                du_y.push_back (nedelec_gradient (neighbor, u, 7));
              break;
            case 6:
              if (node_side == 4)
                du_y.push_back (nedelec_gradient (neighbor, u, 3));
              break;
            case 7:
              if (node_side == 5)
                du_y.push_back (nedelec_gradient (neighbor, u, 7));
              break;
            }

          if (weights_y.size () < du_y.size ())
            {
              weights_y[0] += 2 / hy;
              weights_y.push_back (-1 / hy);
            }
        }
    }

  // If on any boundary/interface normal to z
  if (du_z.size () < 2)
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
          // of the side along z containing "node".
          switch (node)
            {
            case 0:
              node_side = 4;
              break;
            case 1:
              node_side = 5;
              break;
            case 2:
              node_side = 6;
              break;
            case 3:
              node_side = 7;
              break;
            case 4:
              node_side = 0;
              break;
            case 5:
              node_side = 1;
              break;
            case 6:
              node_side = 2;
              break;
            case 7:
              node_side = 3;
              break;
            }

          if (neighbor->gt (node) != quadrant->gt (node_side) ||
              neighbor->is_hanging (node))
            continue;

          hz = neighbor->p (2, 4) - neighbor->p (2, 0);

          switch (node)
            {
            case 0:
              if (node_side == 4)
                du_z.push_back (nedelec_gradient (neighbor, u, 1));
              break;
            case 1:
              if (node_side == 5)
                du_z.push_back (nedelec_gradient (neighbor, u, 4));
              break;
            case 2:
              if (node_side == 6)
                du_z.push_back (nedelec_gradient (neighbor, u, 0));
              break;
            case 3:
              if (node_side == 7)
                du_z.push_back (nedelec_gradient (neighbor, u, 5));
              break;
            case 4:
              if (node_side == 0)
                du_z.push_back (nedelec_gradient (neighbor, u, 1));
              break;
            case 5:
              if (node_side == 1)
                du_z.push_back (nedelec_gradient (neighbor, u, 4));
              break;
            case 6:
              if (node_side == 2)
                du_z.push_back (nedelec_gradient (neighbor, u, 0));
              break;
            case 7:
              if (node_side == 3)
                du_z.push_back (nedelec_gradient (neighbor, u, 5));
              break;
            }

          if (weights_z.size () < du_z.size ())
            {
              weights_z[0] += 2 / hz;
              weights_z.push_back (-1 / hz);
            }
        }
    }

  assert (du_x.size () <= 2 && du_y.size () <= 2
          && du_z.size () <= 2);

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

  if (du_z.size () == 2)
    {
      assigned_z = true;

      for (unsigned int iz = 0; iz < du_z.size (); ++iz)
        du_z_star +=
          du_z[iz] * weights_z[iz];

      du_z_star /=
        std::accumulate (weights_z.begin (),
                         weights_z.end (), 0.0);
    }

  return std::make_tuple (du_x_star, du_y_star, du_z_star,
                          assigned_x, assigned_y, assigned_z);
}

// Specialization ----> std::vector<double>
template <>
gradient3<std::vector<double>>
bim3c_quadtree_pde_recovered_gradient (tmesh_3d& mesh,
                                       const std::vector<double>& u,
                                       active_fun3 is_active)
{
  std::vector<double> du_x_star (mesh.num_global_nodes (), 0);
  std::vector<double> du_y_star (mesh.num_global_nodes (), 0);
  std::vector<double> du_z_star (mesh.num_global_nodes (), 0);

  std::vector<bool> assigned_x (mesh.num_global_nodes (), 0);
  std::vector<bool> assigned_y (mesh.num_global_nodes (), 0);
  std::vector<bool> assigned_z (mesh.num_global_nodes (), 0);

  std::tuple<double, double, double, bool, bool, bool> du_star_loc;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      // Loop over non-hanging vertices of current quadrant.
      for (int node = 0; node < 8; ++node)
        {
          if (quadrant->is_hanging (node) ||
              (assigned_x[quadrant->gt (node)] &&
               assigned_y[quadrant->gt (node)] &&
               assigned_z[quadrant->gt (node)]))
            continue;

          // Skip inactive quadrants.
          if (! is_active(quadrant))
            continue;

          du_star_loc =
            bim3c_recovered_gradient_loc (quadrant, node,
                                          u, is_active);

          du_x_star [quadrant->gt (node)] = std::get<0> (du_star_loc);
          du_y_star [quadrant->gt (node)] = std::get<1> (du_star_loc);
          du_z_star [quadrant->gt (node)] = std::get<2> (du_star_loc);

          assigned_x[quadrant->gt (node)] = std::get<3> (du_star_loc);
          assigned_y[quadrant->gt (node)] = std::get<4> (du_star_loc);
          assigned_z[quadrant->gt (node)] = std::get<5> (du_star_loc);
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

  MPI_Allreduce (MPI_IN_PLACE, du_z_star.data (),
                 du_z_star.size (), MPI_DOUBLE,
                 op, MPI_COMM_WORLD);

  return std::make_tuple (du_x_star, du_y_star, du_z_star);
}

// Specialization ----> distributed_vector
template <>
gradient3<distributed_vector>
bim3c_quadtree_pde_recovered_gradient (tmesh_3d& mesh,
                                       const distributed_vector& u,
                                       active_fun3 is_active)
{
  distributed_vector du_x_star (mesh.num_owned_nodes ());
  distributed_vector du_y_star (mesh.num_owned_nodes ());
  distributed_vector du_z_star (mesh.num_owned_nodes ());

  distributed_vector assigned_x (mesh.num_owned_nodes ());
  distributed_vector assigned_y (mesh.num_owned_nodes ());
  distributed_vector assigned_z (mesh.num_owned_nodes ());

  std::tuple<double, double, double, bool, bool, bool> du_star_loc;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      // Loop over non-hanging vertices of current quadrant.
      for (int node = 0; node < 8; ++node)
        {
          if (assigned_x[quadrant->gt (node)] &&
              assigned_y[quadrant->gt (node)] &&
              assigned_z[quadrant->gt (node)])
            continue;

          // Assemble non-hanging nodes
          if (! quadrant->is_hanging (node))
            {
              // Assemble entries related to inactive quadrants
              du_x_star [quadrant->gt (node)] = 0;
              du_y_star [quadrant->gt (node)] = 0;
              du_z_star [quadrant->gt (node)] = 0;

              // Skip inactive quadrants.
              if (! is_active(quadrant))
                continue;

              du_star_loc =
                bim3c_recovered_gradient_loc (quadrant, node,
                                              u, is_active);

              du_x_star [quadrant->gt (node)] = std::get<0> (du_star_loc);
              du_y_star [quadrant->gt (node)] = std::get<1> (du_star_loc);
              du_z_star [quadrant->gt (node)] = std::get<2> (du_star_loc);

              assigned_x[quadrant->gt (node)] = std::get<3> (du_star_loc);
              assigned_y[quadrant->gt (node)] = std::get<4> (du_star_loc);
              assigned_z[quadrant->gt (node)] = std::get<5> (du_star_loc);
            }
          // Assemble parents
          else
            {
              int np = quadrant->num_parents(node);
              for (int pp = 0; pp < np; ++pp)
                {
                  du_x_star [quadrant->gparent(pp,node)] += 0;
                  du_y_star [quadrant->gparent(pp,node)] += 0;
                  du_z_star [quadrant->gparent(pp,node)] += 0;
                }
            }
        }
    }

  du_x_star.assemble (replace_op);
  du_y_star.assemble (replace_op);
  du_z_star.assemble (replace_op);

  return std::make_tuple (du_x_star, du_y_star, du_z_star);
}

template <class T>
void
compute_solution_if_hanging (std::array<double, 8>& u_star_loc,
                             const gradient3<T>& du,
                             tmesh_3d::quadrant_iterator & quadrant,
                             int n, int pp, int i)
{
  double hx = quadrant->p (0, 7) - quadrant->p (0, 0);
  double hy = quadrant->p (1, 7) - quadrant->p (1, 0);
  double hz = quadrant->p (2, 7) - quadrant->p (2, 0);

  int np = quadrant->num_parents(n);
  int q1 = -1;	// index of parent on the same side
  int q2 = -1;	// index of first parent on the opposite side
  int q3 = -1;  // index of second parent on the opposite side
  if (np > 2)
    {
      if (pp == 0)
        {
          q1 = 1;
          q2 = 2;
          q3 = 3;
        }
      else if (pp == 1)
        {
          q1 = 0;
          q2 = 3;
          q3 = 2;
        }
      else if (pp == 2)
        {
          q1 = 3;
          q2 = 0;
          q3 = 1;
        }
      else if (pp == 3)
        {
          q1 = 2;
          q2 = 1;
          q3 = 0;
        }
    }
  else
    q1 = 1-pp;		// in case of two parents, their indices are always 0,1

  auto normal_to_x =
    [hy, hz, &du, pp, q1, q2, q3, n]
    (tmesh_3d::quadrant_iterator & quadrant, std::array<double, 8>& u_star_loc)
    {
      u_star_loc[n] +=
      hz * (std::get<2>(du)[quadrant->gparent (q1, n)] +
            std::get<2>(du)[quadrant->gparent (pp, n)] -
            std::get<2>(du)[quadrant->gparent (q2, n)] -
            std::get<2>(du)[quadrant->gparent (q3, n)]) / 16;
      u_star_loc[n] +=
      hy * (std::get<1>(du)[quadrant->gparent (pp, n)] +
            std::get<1>(du)[quadrant->gparent (q2, n)] -
            std::get<1>(du)[quadrant->gparent (q1, n)] -
            std::get<1>(du)[quadrant->gparent (q3, n)]) / 16;
    };

  auto normal_to_y =
    [hx, hz, &du, pp, q1, q2, q3, n]
    (tmesh_3d::quadrant_iterator & quadrant, std::array<double, 8>& u_star_loc)
    {
      u_star_loc[n] +=
      hx * (std::get<0>(du)[quadrant->gparent (pp, n)] +
            std::get<0>(du)[quadrant->gparent (q2, n)] -
            std::get<0>(du)[quadrant->gparent (q1, n)] -
            std::get<0>(du)[quadrant->gparent (q3, n)]) / 16;
      u_star_loc[n] +=
      hz * (std::get<2>(du)[quadrant->gparent (pp, n)] +
            std::get<2>(du)[quadrant->gparent (q1, n)] -
            std::get<2>(du)[quadrant->gparent (q2, n)] -
            std::get<2>(du)[quadrant->gparent (q3, n)]) / 16;
    };

  auto normal_to_z =
    [hx, hy, &du, pp, q1, q2, q3, n]
    (tmesh_3d::quadrant_iterator & quadrant, std::array<double, 8>& u_star_loc)
    {
      u_star_loc[n] +=
      hx * (std::get<0>(du)[quadrant->gparent (pp, n)] +
            std::get<0>(du)[quadrant->gparent (q2, n)] -
            std::get<0>(du)[quadrant->gparent (q1, n)] -
            std::get<0>(du)[quadrant->gparent (q3, n)]) / 16;
      u_star_loc[n] +=
      hy * (std::get<1>(du)[quadrant->gparent (pp, n)] +
            std::get<1>(du)[quadrant->gparent (q1, n)] -
            std::get<1>(du)[quadrant->gparent (q2, n)] -
            std::get<1>(du)[quadrant->gparent (q3, n)]) / 16;
    };

  if (n == 0)	// ------- node 0
    {
      if (np > 2) // four parents
        {
          if (i == 2) // face normal to x
            normal_to_x(quadrant,u_star_loc);
          else if (i == 8)  // face normal to z
            normal_to_z(quadrant,u_star_loc);
          else if (i == 9)  // face normal to y
            normal_to_y(quadrant,u_star_loc);
        }
      else  // two parents
        {
          if (i == 1) // edge directed along x
            u_star_loc[n] +=
              (2 * hx) * (std::get<0>(du)[quadrant->gparent (q1, n)] -
                          std::get<0>(du)[quadrant->gparent (pp, n)]) / 8;
          else if (i == 2) // edge directed along y
            u_star_loc[n] +=
              (2 * hy) * (std::get<1>(du)[quadrant->gparent (q1, n)] -
                          std::get<1>(du)[quadrant->gparent (pp, n)]) / 8;
          else if (i == 4)  // edge directed along z
            u_star_loc[n] +=
              (2 * hz) * (std::get<2>(du)[quadrant->gparent (q1, n)] -
                          std::get<2>(du)[quadrant->gparent (pp, n)]) / 8;
        }
    }
  else if (n == 1)	// ------- node 1
    {
      if (np > 2) // four parents
        {
          if (i == 3) // face normal to x
            normal_to_x(quadrant,u_star_loc);
          else if (i == 8)  // face normal to z
            normal_to_z(quadrant,u_star_loc);
          else if (i == 9)  // face normal to y
            normal_to_y(quadrant,u_star_loc);
        }
      else  // two parents
        {
          if (i == 0) // edge directed along x
            u_star_loc[n] +=
              (2 * hx) * (std::get<0>(du)[quadrant->gparent (pp, n)] -
                          std::get<0>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 3)  // edge directed along y
            u_star_loc[n] +=
              (2 * hy) * (std::get<1>(du)[quadrant->gparent (q1, n)] -
                          std::get<1>(du)[quadrant->gparent (pp, n)]) / 8;
          else if (i == 5)  // edge directed along z
            u_star_loc[n] +=
              (2 * hz) * (std::get<2>(du)[quadrant->gparent (q1, n)] -
                          std::get<2>(du)[quadrant->gparent (pp, n)]) / 8;
        }
    }
  else if (n == 2)	// ------- node 2
    {
      if (np > 2) // four parents
        {
          if (i == 3) // face normal to y
            normal_to_y(quadrant,u_star_loc);
          else if (i == 8)  // face normal to z
            normal_to_z(quadrant,u_star_loc);
          else if (i == 9)  // face normal to x
            normal_to_x(quadrant,u_star_loc);
        }
      else  // two parents
        {
          if (i == 0) // edge directed along y
            u_star_loc[n] +=
              (2 * hy) * (std::get<1>(du)[quadrant->gparent (pp, n)] -
                          std::get<1>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 3)  // edge directed along x
            u_star_loc[n] +=
              (2 * hx) * (std::get<0>(du)[quadrant->gparent (q1, n)] -
                          std::get<0>(du)[quadrant->gparent (pp, n)]) / 8;
          else if (i == 6)  // edge directed along z
            u_star_loc[n] +=
              (2 * hz) * (std::get<2>(du)[quadrant->gparent (q1, n)] -
                          std::get<2>(du)[quadrant->gparent (pp, n)]) / 8;
        }
    }
  else if (n == 3)	// ------- node 2
    {
      if (np > 2) // four parents
        {
          if (i == 0) // face normal to z
            normal_to_z(quadrant,u_star_loc);
          else if (i == 1)  // face normal to x
            normal_to_x(quadrant,u_star_loc);
          else if (i == 2)  // face normal to y
            normal_to_y(quadrant,u_star_loc);
        }
      else  // two parents
        {
          if (i == 1) // edge directed along y
            u_star_loc[n] +=
              (2 * hy) * (std::get<1>(du)[quadrant->gparent (pp, n)] -
                          std::get<1>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 2)  // edge directed along x
            u_star_loc[n] +=
              (2 * hx) * (std::get<0>(du)[quadrant->gparent (pp, n)] -
                          std::get<0>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 7)  // edge directed along z
            u_star_loc[n] +=
              (2 * hz) * (std::get<2>(du)[quadrant->gparent (q1, n)] -
                          std::get<2>(du)[quadrant->gparent (pp, n)]) / 8;
        }
    }
  else if (n == 4)	// ------- node 4
    {
      if (np > 2) // four parents
        {
          if (i == 5) // face normal to z
            normal_to_z(quadrant,u_star_loc);
          else if (i == 8)  // face normal to y
            normal_to_y(quadrant,u_star_loc);
          else if (i == 9)  // face normal to x
            normal_to_x(quadrant,u_star_loc);
        }
      else  // two parents
        {
          if (i == 0) // edge directed along z
            u_star_loc[n] +=
              (2 * hz) * (std::get<2>(du)[quadrant->gparent (pp, n)] -
                          std::get<2>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 5)  // edge directed along x
            u_star_loc[n] +=
              (2 * hx) * (std::get<0>(du)[quadrant->gparent (q1, n)] -
                          std::get<0>(du)[quadrant->gparent (pp, n)]) / 8;
          else if (i == 6)  // edge directed along y
            u_star_loc[n] +=
              (2 * hy) * (std::get<1>(du)[quadrant->gparent (q1, n)] -
                          std::get<1>(du)[quadrant->gparent (pp, n)]) / 8;
        }
    }
  else if (n == 5)	// ------- node 5
    {
      if (np > 2) // four parents
        {
          if (i == 0)	// face normal to y
            normal_to_y(quadrant,u_star_loc);
          else if (i == 1)  // face normal to x
            normal_to_x(quadrant,u_star_loc);
          else if (i == 4)  // face normal to z
            normal_to_z(quadrant,u_star_loc);
        }
      else  // two parents
        {
          if (i == 1) // edge directed along z
            u_star_loc[n] +=
              (2 * hz) * (std::get<2>(du)[quadrant->gparent (pp, n)] -
                          std::get<2>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 4)  // edge directed along x
            u_star_loc[n] +=
              (2 * hx) * (std::get<0>(du)[quadrant->gparent (pp, n)] -
                          std::get<0>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 7)  // edge directed along y
            u_star_loc[n] +=
              (2 * hy) * (std::get<1>(du)[quadrant->gparent (q1, n)] -
                          std::get<1>(du)[quadrant->gparent (pp, n)]) / 8;
        }
    }
  else if (n == 6)	// ------- node 6
    {
      if (np > 2) // four parents
        {
          if (i == 2) // face normal to y
            normal_to_y(quadrant,u_star_loc);
          else if (i == 4)  // face normal to z
            normal_to_z(quadrant,u_star_loc);
        }
      else  // two parents
        {
          if (i == 2) // edge directed along z
            u_star_loc[n] +=
              (2 * hz) * (std::get<2>(du)[quadrant->gparent (pp, n)] -
                          std::get<2>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 4)  // edge directed along y
            u_star_loc[n] +=
              (2 * hy) * (std::get<1>(du)[quadrant->gparent (pp, n)] -
                          std::get<1>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 7)  // edge directed along x
            u_star_loc[n] +=
              (2 * hx) * (std::get<0>(du)[quadrant->gparent (q1, n)] -
                          std::get<0>(du)[quadrant->gparent (pp, n)]) / 8;
        }
    }
  else if (n == 7)	// ------- node 7
    {
      if (np > 2)	// four parents;
        {
          if (i == 1)	// face normal to x
            normal_to_x(quadrant,u_star_loc);
          else if (i == 2)	// face normal to y
            normal_to_y(quadrant,u_star_loc);
          else if (i == 4)	// face normal to z
            normal_to_z(quadrant,u_star_loc);
        }
      else	// two parents;
        {
          if (i == 3)	// edge directed along z
            u_star_loc[n] +=
              (2 * hz) * (std::get<2>(du)[quadrant->gparent (pp, n)] -
                          std::get<2>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 5)	// edge directed along y
            u_star_loc[n] +=
              (2 * hy) * (std::get<1>(du)[quadrant->gparent (pp, n)] -
                          std::get<1>(du)[quadrant->gparent (q1, n)]) / 8;
          else if (i == 6)	// edge directed along x
            u_star_loc[n] +=
              (2 * hx) * (std::get<0>(du)[quadrant->gparent (pp, n)] -
                          std::get<0>(du)[quadrant->gparent (q1, n)]) / 8;
        }
    }
}

template <class T>
q2_vec3
bim3c_quadtree_pde_recovered_solution (tmesh_3d& mesh,
                                       const T& u,
                                       const gradient3<T>& du)
{
  q2_vec3 u_star (mesh.num_local_quadrants (),
                  std::array<double, 27> ({0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0}));

  double hx = 0.0, hy = 0.0, hz = 0.0;

  std::array<double, 8> u_star_loc,
    du_x_star_loc,
    du_y_star_loc,
    du_z_star_loc;


  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p (0, 7) - quadrant->p (0, 0);
      hy = quadrant->p (1, 7) - quadrant->p (1, 0);
      hz = quadrant->p (2, 7) - quadrant->p (2, 0);

      // Compute values at vertices.
      for (int n = 0; n < 8; ++n)
        {
          if (! quadrant->is_hanging (n))
            {
              u_star_loc[n] = u[quadrant->gt (n)];
              du_x_star_loc[n] = std::get<0>(du)[quadrant->gt (n)];
              du_y_star_loc[n] = std::get<1>(du)[quadrant->gt (n)];
              du_z_star_loc[n] = std::get<2>(du)[quadrant->gt (n)];
            }
          else
            {
              // value of u_star_loc, du_x_star_loc, du_y_star_loc, and
              // du_z_star_loc are computed as the average over all the
              // parents
              int np = quadrant->num_parents(n);
              u_star_loc[n] = 0.;
              du_x_star_loc[n] = 0.;
              du_y_star_loc[n] = 0.;
              du_z_star_loc[n] = 0.;
              for (int pp = 0; pp < np; ++pp)
                {
                  u_star_loc[n] += u[quadrant->gparent(pp,n)];
                  du_x_star_loc[n] += std::get<0>(du)
                    [quadrant->gparent (pp, n)];
                  du_y_star_loc[n] += std::get<1>(du)
                    [quadrant->gparent (pp, n)];
                  du_z_star_loc[n] += std::get<2>(du)
                    [quadrant->gparent (pp, n)];
                }
              u_star_loc[n] /= np;
              du_x_star_loc[n] /= np;
              du_y_star_loc[n] /= np;
              du_z_star_loc[n] /= np;

              // Determine whether n is hanging on an edge
              // directed along the x or y or z direction.
              int i = 0; int p = 0;
              bool found = false;

              // we loop over the np parents of node n, looking for the index
              // of the first parent whose global index is equal to the index
              // of the current node
              for (; i < 8 && !found; ++i)
                {
                  for (int pp = 0; pp < np && !found; ++pp)
                    {
                      int par = quadrant->parent(pp,n);
                      int tt = quadrant->t(i);
                      if (i != n && par == tt)
                        {
                          if (np > 2)
                            {
                              if (i == 0)
                                {
                                  found = false;
                                  int j = 1;
                                  for (; j < 8 && !found; ++j)
                                    {
                                      for (int pp = 0; pp < np && !found; ++pp)
                                        {
                                          int par = quadrant->parent(pp,n);
                                          int tt = quadrant->t(j);
                                          if (j!=n && par == tt)
                                            {
                                              found = true;
                                              --j;
                                            }
                                        }
                                    }
                                  if (n == 1)
                                    {
                                      // node 1 --->
                                      // i=0,j=2 -> face normal to z
                                      // i=0,j=4 -> face normal to y
                                      if (j == 2)
                                        i = 9;
                                      else if (j == 4)
                                        i = 10;
                                    }
                                  else if (n == 2)
                                    {
                                      // node 2 --->
                                      // i=0,j=1 -> face normal to z
                                      // i=0,j=4 -> face normal to x
                                      if (j == 1)
                                        i = 9;
                                      else if (j == 4)
                                        i = 10;
                                    }
                                  else if (n == 4)
                                    {
                                      // node 4 --->
                                      // i=0,j=1 -> face normal to y
                                      // i=0,j=2 -> face normal to x
                                      if (j == 1)
                                        i = 9;
                                      else if (j == 2)
                                        i = 10;
                                    }
                                }
                              else if (n == 0 && i == 1)
                                {
                                  // node 0 --->
                                  // i=1,j=2 -> face normal to z
                                  // i=1,j=4 -> face normal to y
                                  found = false;
                                  int j = 2;
                                  for (; j < 8 && !found; ++j)
                                    {
                                      for (int pp = 0; pp < np && !found; ++pp)
                                        {
                                          int par = quadrant->parent(pp,n);
                                          int tt = quadrant->t(j);
                                          if (j!=n && par == tt)
                                            {
                                              found = true;
                                              --j;
                                            }
                                        }
                                    }
                                  if (j == 2)
                                    i = 9;
                                  else if (j == 4)
                                    i = 10;
                                }
                            }
                          p = pp;
                          found = true;
                          --i;
                        }
                    }
                }

              // Compute recovered solution at the
              // double-sized neighbor element.
              compute_solution_if_hanging (u_star_loc,du,quadrant,n,p,i);
            }

          u_star[quadrant->get_forest_quad_idx ()][n] = u_star_loc[n];
        }

      // Compute values at edges
      // face 0
      u_star[quadrant->get_forest_quad_idx ()][8] =
        0.5 * (u_star_loc[2] + u_star_loc[6])
        + hz * (du_z_star_loc[2] - du_z_star_loc[6]) / 8;

      u_star[quadrant->get_forest_quad_idx ()][9] =
        0.5 * (u_star_loc[0] + u_star_loc[4])
        + hz * (du_z_star_loc[0] - du_z_star_loc[4]) / 8;

      u_star[quadrant->get_forest_quad_idx ()][10] =
        0.5 * (u_star_loc[2] + u_star_loc[0])
        + hy * (du_y_star_loc[0] - du_y_star_loc[2]) / 8;

      u_star[quadrant->get_forest_quad_idx ()][11] =
        0.5 * (u_star_loc[6] + u_star_loc[4])
        + hy * (du_y_star_loc[4] - du_y_star_loc[6]) / 8;

      // face 1
      u_star[quadrant->get_forest_quad_idx ()][12] =
        0.5 * (u_star_loc[1] + u_star_loc[5])
        + hz * (du_z_star_loc[1] - du_z_star_loc[5]) / 8;

      u_star[quadrant->get_forest_quad_idx ()][13] =
        0.5 * (u_star_loc[3] + u_star_loc[7])
        + hz * (du_z_star_loc[3] - du_z_star_loc[7]) / 8;

      u_star[quadrant->get_forest_quad_idx ()][14] =
        0.5 * (u_star_loc[1] + u_star_loc[3])
        + hy * (du_y_star_loc[1] - du_y_star_loc[3]) / 8;

      u_star[quadrant->get_forest_quad_idx ()][15] =
        0.5 * (u_star_loc[5] + u_star_loc[7])
        + hy * (du_y_star_loc[5] - du_y_star_loc[7]) / 8;

      // face 2
      u_star[quadrant->get_forest_quad_idx ()][16] =
        0.5 * (u_star_loc[0] + u_star_loc[1])
        + hx * (du_x_star_loc[0] - du_x_star_loc[1]) / 8;

      u_star[quadrant->get_forest_quad_idx ()][17] =
        0.5 * (u_star_loc[4] + u_star_loc[5])
        + hx * (du_x_star_loc[4] - du_x_star_loc[5]) / 8;

      // face 3
      u_star[quadrant->get_forest_quad_idx ()][18] =
        0.5 * (u_star_loc[3] + u_star_loc[2])
        + hx * (du_x_star_loc[2] - du_x_star_loc[3]) / 8;

      u_star[quadrant->get_forest_quad_idx ()][19] =
        0.5 * (u_star_loc[7] + u_star_loc[6])
        + hx * (du_x_star_loc[6] - du_x_star_loc[7]) / 8;

      // Compute values at faces midpoints.
      // face 0
      u_star[quadrant->get_forest_quad_idx ()][20] =
        0.25 * (u_star[quadrant->get_forest_quad_idx ()][8] +
                u_star[quadrant->get_forest_quad_idx ()][9] +
                u_star[quadrant->get_forest_quad_idx ()][10] +
                u_star[quadrant->get_forest_quad_idx ()][11])
        + hz * 0.5 * (du_z_star_loc[2] + du_z_star_loc[0] -
                      du_z_star_loc[6] - du_z_star_loc[4]) / 16
        + hy * 0.5 * (du_y_star_loc[0] + du_y_star_loc[4] -
                      du_y_star_loc[2] - du_y_star_loc[6]) / 16;

      // face 1
      u_star[quadrant->get_forest_quad_idx ()][21] =
        0.25 * (u_star[quadrant->get_forest_quad_idx ()][12] +
                u_star[quadrant->get_forest_quad_idx ()][13] +
                u_star[quadrant->get_forest_quad_idx ()][14] +
                u_star[quadrant->get_forest_quad_idx ()][15])
        + hz * 0.5 * (du_z_star_loc[1] + du_z_star_loc[3] -
                      du_z_star_loc[5] - du_z_star_loc[7]) / 16
        + hy * 0.5 * (du_y_star_loc[1] + du_y_star_loc[5] -
                      du_y_star_loc[3] - du_y_star_loc[7]) / 16;

      // face 2
      u_star[quadrant->get_forest_quad_idx ()][22] =
        0.25 * (u_star[quadrant->get_forest_quad_idx ()][9] +
                u_star[quadrant->get_forest_quad_idx ()][12] +
                u_star[quadrant->get_forest_quad_idx ()][16] +
                u_star[quadrant->get_forest_quad_idx ()][17])
        + hx * 0.5 * (du_x_star_loc[0] + du_x_star_loc[4] -
                      du_x_star_loc[1] - du_x_star_loc[5]) / 16
        + hz * 0.5 * (du_z_star_loc[0] + du_z_star_loc[1] -
                      du_z_star_loc[4] - du_z_star_loc[5]) / 16;

      // face 3
      u_star[quadrant->get_forest_quad_idx ()][23] =
        0.25 * (u_star[quadrant->get_forest_quad_idx ()][8] +
                u_star[quadrant->get_forest_quad_idx ()][13] +
                u_star[quadrant->get_forest_quad_idx ()][18] +
                u_star[quadrant->get_forest_quad_idx ()][19])
        + hx * 0.5 * (du_x_star_loc[2] + du_x_star_loc[6] -
                      du_x_star_loc[3] - du_x_star_loc[7]) / 16
        + hz * 0.5 * (du_z_star_loc[3] + du_z_star_loc[2] -
                      du_z_star_loc[6] - du_z_star_loc[7]) / 16;

      // face 4
      u_star[quadrant->get_forest_quad_idx ()][24] =
        0.25 * (u_star[quadrant->get_forest_quad_idx ()][10] +
                u_star[quadrant->get_forest_quad_idx ()][14] +
                u_star[quadrant->get_forest_quad_idx ()][16] +
                u_star[quadrant->get_forest_quad_idx ()][18])
        + hx * 0.5 * (du_x_star_loc[0] + du_x_star_loc[2] -
                      du_x_star_loc[1] - du_x_star_loc[3]) / 16
        + hy * 0.5 * (du_y_star_loc[0] + du_y_star_loc[1] -
                      du_y_star_loc[2] - du_y_star_loc[3]) / 16;

      // face 5
      u_star[quadrant->get_forest_quad_idx ()][25] =
        0.25 * (u_star[quadrant->get_forest_quad_idx ()][11] +
                u_star[quadrant->get_forest_quad_idx ()][15] +
                u_star[quadrant->get_forest_quad_idx ()][17] +
                u_star[quadrant->get_forest_quad_idx ()][19])
        + hx * 0.5 * (du_x_star_loc[4] + du_x_star_loc[6] -
                      du_x_star_loc[5] - du_x_star_loc[7]) / 16
        + hy * 0.5 * (du_y_star_loc[4] + du_y_star_loc[5] -
                      du_y_star_loc[6] - du_y_star_loc[7]) / 16;

      // Compute value at cell midpoint.
      u_star[quadrant->get_forest_quad_idx ()][26] =
        (u_star[quadrant->get_forest_quad_idx ()][20] +
         u_star[quadrant->get_forest_quad_idx ()][21] +
         u_star[quadrant->get_forest_quad_idx ()][22] +
         u_star[quadrant->get_forest_quad_idx ()][23] +
         u_star[quadrant->get_forest_quad_idx ()][24] +
         u_star[quadrant->get_forest_quad_idx ()][25]) / 6
        + hx * 0.25 * (du_x_star_loc[0] + du_x_star_loc[2] +
                       du_x_star_loc[4] + du_x_star_loc[6] -
                       du_x_star_loc[1] - du_x_star_loc[3] -
                       du_x_star_loc[5] - du_x_star_loc[7]) / 24
        + hy * 0.25 * (du_y_star_loc[0] + du_y_star_loc[1] +
                       du_y_star_loc[4] + du_y_star_loc[5] -
                       du_y_star_loc[2] - du_y_star_loc[3] -
                       du_y_star_loc[6] - du_y_star_loc[7]) / 24
        + hz * 0.25 * (du_z_star_loc[0] + du_z_star_loc[1] +
                       du_z_star_loc[2] + du_z_star_loc[3] -
                       du_z_star_loc[4] - du_z_star_loc[5] -
                       du_z_star_loc[6] - du_z_star_loc[7]) / 24;
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


// Approximate integral of fun on [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]].
static double
quad_integral (const double *x, const double *y, const double *z,
               std::function<double (double, double, double)> fun)
{
  int ix, jy, kz;
  double wx = 0, nx = 0, sum = 0;
  for (ix = 0; ix < 4; ++ix)
    {
      wx = xformw (x, gw[ix]);
      nx = xformx (x, gn[ix]);
      double wy = 0, ny = 0;
      for (jy = 0; jy < 4; ++jy)
        {
          wy = xformw (y, gw[jy]);
          ny = xformx (y, gn[jy]);
          for (kz = 0; kz < 4; ++kz)
            sum += fun (nx, ny, xformx(z, gn[kz])) *
              wx * wy * xformw(z, gw[kz]);
        }
    }
  return (sum);
}

// Evaluate Nedelec x-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
dudx (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double db1 = (u[1] - u[0]);
  double dt1 = (u[3] - u[2]);

  double db2 = (u[5] - u[4]);
  double dt2 = (u[7] - u[6]);

  double result = (db1 * (y[1] - Y) * (z[1] - Z) +
                   dt1 * (Y - y[0]) * (z[1] - Z) +
                   db2 * (y[1] - Y) * (Z - z[0]) +
                   dt2 * (Y - y[0]) * (Z - z[0])) / hxhyhz;

  return  result;
}

// Evaluate Nedelec y-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
dudy (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double dl1 = (u[2] - u[0]);
  double dr1 = (u[3] - u[1]);

  double dl2 = (u[6] - u[4]);
  double dr2 = (u[7] - u[5]);

  double result = (dl1 * (x[1] - X) * (z[1] - Z) +
                   dr1 * (X - x[0]) * (z[1] - Z) +
                   dl2 * (x[1] - X) * (Z - z[0]) +
                   dr2 * (X - x[0]) * (Z - z[0])) / hxhyhz;

  return result;
}

// Evaluate Nedelec z-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
dudz (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double d11 = (u[4] - u[0]);
  double d21 = (u[5] - u[1]);

  double d12 = (u[6] - u[2]);
  double d22 = (u[7] - u[3]);

  double result = (d11 * (x[1] - X) * (y[1] - Y) +
                   d21 * (X - x[0]) * (y[1] - Y) +
                   d12 * (x[1] - X) * (Y - y[0]) +
                   d22 * (X - x[0]) * (Y - y[0])) / hxhyhz;

  return result;
}

// Evaluate u (using Q1 basis functions
// on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
q1 (double X, double Y, double Z, const double *x,
    const double *y, const double *z, const double *u)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);
  double Xx1 = (X - x[1]);
  double Xx0 = (X - x[0]);
  double Yy1 = (Y - y[1]);
  double Yy0 = (Y - y[0]);
  double Zz1 = (Z - z[1]);
  double Zz0 = (Z - z[0]);

  double num = u[0] * -Xx1  * Yy1 * Zz1 +
    u[1] * Xx0 * Yy1 * Zz1 +
    u[2] * Xx1 * Yy0 * Zz1 +
    u[3] * -Xx0  * Yy0 * Zz1 +
    u[4] * Xx1  * Yy1 * Zz0 +
    u[5] * -Xx0 * Yy1 * Zz0 +
    u[6] * -Xx1 * Yy0 * Zz0 +
    u[7] * Xx0  * Yy0 * Zz0;

  return (num / hxhyhz);
}

// Evaluate u (using Q2 basis functions
// on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
q2 (double X, double Y, double Z, const double *x,
    const double *y, const double *z, const double *u)
{
  double xc = 0.5 * (x[0] + x[1]);
  double yc = 0.5 * (y[0] + y[1]);
  double zc = 0.5 * (z[0] + z[1]);

  double Xx1 = (X - x[1]);
  double Xx0 = (X - x[0]);
  double Xxc = (X - xc);
  double Yy1 = (Y - y[1]);
  double Yy0 = (Y - y[0]);
  double Yyc = (Y - yc);
  double Zz1 = (Z - z[1]);
  double Zz0 = (Z - z[0]);
  double Zzc = (Z - zc);

  double hxhyhz2 = std::pow ((x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]), 2);

  double num = (u[0] * 8  * Xxc * Xx1 * Yyc * Yy1 * Zzc * Zz1 +
                u[1] * 8  * Xx0 * Xxc * Yyc * Yy1 * Zzc * Zz1 +
                u[2] * 8  * Xxc * Xx1 * Yy0 * Yyc * Zzc * Zz1 +
                u[3] * 8  * Xx0 * Xxc * Yy0 * Yyc * Zzc * Zz1 +

                u[4] * 8  * Xxc * Xx1 * Yyc * Yy1 * Zz0 * Zzc +
                u[5] * 8  * Xx0 * Xxc * Yyc * Yy1 * Zz0 * Zzc +
                u[6] * 8  * Xxc * Xx1 * Yy0 * Yyc * Zz0 * Zzc +
                u[7] * 8  * Xx0 * Xxc * Yy0 * Yyc * Zz0 * Zzc +

                u[8] * -16  * Xxc * Xx1 * Yy0 * Yyc * Zz0 * Zz1 +
                u[9] * -16  * Xxc * Xx1 * Yyc * Yy1 * Zz0 * Zz1 +
                u[10] * -16  * Xxc * Xx1 * Yy0 * Yy1 * Zzc * Zz1 +
                u[11] * -16  * Xxc * Xx1 * Yy0 * Yy1 * Zz0 * Zzc +

                u[12] * -16  * Xx0 * Xxc * Yyc * Yy1 * Zz0 * Zz1 +
                u[13] * -16  * Xx0 * Xxc * Yy0 * Yyc * Zz0 * Zz1 +
                u[14] * -16  * Xx0 * Xxc * Yy0 * Yy1 * Zzc * Zz1 +
                u[15] * -16  * Xx0 * Xxc * Yy0 * Yy1 * Zz0 * Zzc +

                u[16] * -16  * Xx0 * Xx1 * Yyc * Yy1 * Zzc * Zz1 +
                u[17] * -16  * Xx0 * Xx1 * Yyc * Yy1 * Zz0 * Zzc +

                u[18] * -16  * Xx0 * Xx1 * Yy0 * Yyc * Zzc * Zz1 +
                u[19] * -16  * Xx0 * Xx1 * Yy0 * Yyc * Zz0 * Zzc +

                u[20] * 32  * Xxc * Xx1 * Yy0 * Yy1 * Zz0 * Zz1 +

                u[21] * 32  * Xx0 * Xxc * Yy0 * Yy1 * Zz0 * Zz1 +

                u[22] * 32  * Xx0 * Xx1 * Yyc * Yy1 * Zz0 * Zz1 +

                u[23] * 32  * Xx0 * Xx1 * Yy0 * Yyc * Zz0 * Zz1 +

                u[24] * 32  * Xx0 * Xx1 * Yy0 * Yy1 * Zzc * Zz1 +

                u[25] * 32  * Xx0 * Xx1 * Yy0 * Yy1 * Zz0 * Zzc +

                u[26] * -64 * Xx0 * Xx1 * Yy0 * Yy1 * Zz0 * Zz1);

  return num / hxhyhz2;
}

// Compute ||grad^* u - grad u||_L^2(q).
template <class T>
double
estimator_grad (tmesh_3d::quadrant_iterator q,
                const gradient3<T>& du_star,
                const T& u)
{
  double
    x[2] = {q->p(0,0), q->p(0,7)},
    y[2] = {q->p(1,0), q->p(1,7)},
    z[2] = {q->p(2,0), q->p(2,7)};

  double dudxstar_loc[8] = {0,0,0,0,0,0,0,0};
  double dudystar_loc[8] = {0,0,0,0,0,0,0,0};
  double dudzstar_loc[8] = {0,0,0,0,0,0,0,0};

  double u_loc[8] = {0,0,0,0,0,0,0,0};

  for (int ii = 0; ii < 8; ++ii)
    if (! q->is_hanging (ii))
      {
        dudxstar_loc[ii] = (std::get<0>(du_star))[q->gt (ii)];
        dudystar_loc[ii] = (std::get<1>(du_star))[q->gt (ii)];
        dudzstar_loc[ii] = (std::get<2>(du_star))[q->gt (ii)];
        u_loc[ii] = u[q->gt(ii)];
      }
    else
      {
        int np = q->num_parents(ii);
        for (int pp = 0; pp < np; ++pp)
          {
            dudxstar_loc[ii] += (std::get<0>(du_star))[q->gparent (pp, ii)];
            dudystar_loc[ii] += (std::get<1>(du_star))[q->gparent (pp, ii)];
            dudzstar_loc[ii] += (std::get<2>(du_star))[q->gparent (pp, ii)];
            u_loc[ii] += u[q->gparent (pp, ii)];
          }
        dudxstar_loc[ii] /= np;
        dudystar_loc[ii] /= np;
        dudzstar_loc[ii] /= np;
        u_loc[ii] /= np;
      }

  auto fun =
    [&x, &y, &z, &dudxstar_loc, &dudystar_loc, &dudzstar_loc, &u_loc]
    (double X, double Y, double Z) -> double
    {
      return
      std::pow (dudx (X, Y, Z, x, y, z, u_loc) -
                q1 (X, Y, Z, x, y, z, dudxstar_loc), 2) +
      std::pow (dudy (X, Y, Z, x, y, z, u_loc) -
                q1 (X, Y, Z, x, y, z, dudystar_loc), 2) +
      std::pow (dudz (X, Y, Z, x, y, z, u_loc) -
                q1 (X, Y, Z, x, y, z, dudzstar_loc), 2);
    };

  return std::sqrt (quad_integral (x, y, z, fun));
}


template <class T>
int
zz_marker_grad (tmesh_3d::quadrant_iterator q,
                const gradient3<T>& du_star,
                const T& u,
                double limit)
{
  return estimator_grad (q, du_star, u) > limit ? 1 : 0;
}


// Compute ||u^* - u||_L^2(q).
template <class T>
double
estimator_sol (tmesh_3d::quadrant_iterator q,
               const q2_vec3 & ustar,
               const T & u)
{
  double
    x[2] = {q->p(0,0), q->p (0,7)},
    y[2] = {q->p(1,0), q->p (1,7)},
    z[2] = {q->p(2,0), q->p (2,7)};

  double ustar_loc[27] = {0,0,0,0,0,0,0,0,0,
                          0,0,0,0,0,0,0,0,0,
                          0,0,0,0,0,0,0,0,0};
  double u_loc[8] = {0,0,0,0,0,0,0,0};

  for (int ii = 0; ii < 27; ++ii)
    ustar_loc[ii] = (ustar[q->get_forest_quad_idx ()])[ii];

  for (int ii = 0; ii < 8; ++ii)
    {
      if (! q->is_hanging (ii))
        u_loc[ii] = u[q->gt (ii)];
      else
        {
          int np = q->num_parents(ii);
          for (int pp = 0; pp < np; ++pp)
            u_loc[ii] += u[q->gparent (pp, ii)];
          u_loc[ii] /= np;
        }
    }

  auto fun =
    [&x, &y, &z, &ustar_loc, &u_loc]
    (double X, double Y, double Z) -> double
    {
      return
      std::pow (q1 (X, Y, Z, x, y, z, u_loc) -
                q2 (X, Y, Z, x, y, z, ustar_loc), 2);
    };

  return std::sqrt (quad_integral (x, y, z, fun));
}


template <class T>
int
zz_marker_sol (tmesh_3d::quadrant_iterator q,
               const q2_vec3& ustar,
               const T& u,
               double limit)
{
  return estimator_sol (q, ustar, u) > limit ? 1 : 0;
}


// Compute ||u - u_ex||_L^2(q).
template <class T>
double
l2_error (tmesh_3d::quadrant_iterator q,
          const func3 & u_ex,
          const T & u)
{
  double
    x[2] = {q->p (0,0), q->p (0,7)},
    y[2] = {q->p (1,0), q->p (1,7)},
    z[2] = {q->p (2,0), q->p (2,7)};

  double u_loc[8] = {0,0,0,0,0,0,0,0};

  for (int ii = 0; ii < 8; ++ii)
    if (! q->is_hanging (ii))
      u_loc[ii] = u[q->gt (ii)];
    else
      {
        int np = q->num_parents(ii);
        for (int pp = 0; pp < np; ++pp)
          u_loc[ii] += u[q->gparent (pp, ii)];
        u_loc[ii] /= np;
      }

  auto fun =
    [&x, &y, &z, &u_loc, &u_ex]
    (double X, double Y, double Z) -> double
    { return std::pow (q1 (X, Y, Z, x, y, z, u_loc) - u_ex (X, Y, Z), 2); };

  return std::sqrt (quad_integral (x, y, z, fun));
}


// Compute |u - u_ex|_H^1(q).
template <class T>
double
semih1_error (tmesh_3d::quadrant_iterator q,
              const func3 & dudx_ex,
              const func3 & dudy_ex,
              const func3 & dudz_ex,
              const T & u)
{
  double
    x[2] = {q->p (0,0), q->p (0,7)},
    y[2] = {q->p (1,0), q->p (1,7)},
    z[2] = {q->p (2,0), q->p (2,7)};

  double u_loc[8] = {0,0,0,0,0,0,0,0};

  for (int ii = 0; ii < 8; ++ii)
    if (! q->is_hanging (ii))
      u_loc[ii] = u[q->gt (ii)];
    else
      {
        int np = q->num_parents(ii);
        for (int pp = 0; pp < np; ++pp)
          u_loc[ii] += u[q->gparent (pp, ii)];
        u_loc[ii] /= np;
      }

  auto fun =
    [&x, &y, &z, &dudx_ex, &dudy_ex, &dudz_ex, &u_loc]
    (double X, double Y, double Z) -> double
    {
      return
      std::pow (dudx (X, Y, Z, x, y, z, u_loc) -
                dudx_ex (X, Y, Z), 2) +
      std::pow (dudy (X, Y, Z, x, y, z, u_loc) -
                dudy_ex (X, Y, Z), 2) +
      std::pow (dudz (X, Y, Z, x, y, z, u_loc) -
                dudz_ex (X, Y, Z), 2);
    };

  return std::sqrt (quad_integral (x, y, z, fun));
}

// Compute ||u_star - u_ex||_L^2(q).
double
l2_star_error (tmesh_3d::quadrant_iterator q,
               const func3 & u_ex,
               const q2_vec3 & ustar)
{
  double
    x[2] = {q->p (0,0), q->p (0,7)},
    y[2] = {q->p (1,0), q->p (1,7)},
    z[2] = {q->p (2,0), q->p (2,7)};

  double ustar_loc[27] = {0,0,0,0,0,0,0,0,0,
                          0,0,0,0,0,0,0,0,0,
                          0,0,0,0,0,0,0,0,0};

  for (int ii = 0; ii < 27; ++ii)
    ustar_loc[ii] = (ustar[q->get_forest_quad_idx ()])[ii];

  auto fun =
    [&x, &y, &z, &ustar_loc, &u_ex]
    (double X, double Y, double Z) -> double
    {
      return
      std::pow (q2 (X, Y, Z, x, y, z, ustar_loc) - u_ex (X, Y, Z), 2);
    };

  return std::sqrt (quad_integral (x, y, z, fun));
}


// Compute ||du_star - grad(u_ex)||_L^2(q).
template <class T>
double
semih1_star_error (tmesh_3d::quadrant_iterator q,
                   const func3 & dudx_ex,
                   const func3 & dudy_ex,
                   const func3 & dudz_ex,
                   const gradient3<T> & du_star)
{
  double
    x[2] = {q->p (0,0), q->p (0,7)},
    y[2] = {q->p (1,0), q->p (1,7)},
    z[2] = {q->p (2,0), q->p (2,7)};

  double dudxstar_loc[8] = {0,0,0,0,0,0,0,0};
  double dudystar_loc[8] = {0,0,0,0,0,0,0,0};
  double dudzstar_loc[8] = {0,0,0,0,0,0,0,0};

  for (int ii = 0; ii < 8; ++ii)
    if (! q->is_hanging (ii))
      {
        dudxstar_loc[ii] = (std::get<0>(du_star))[q->gt (ii)];
        dudystar_loc[ii] = (std::get<1>(du_star))[q->gt (ii)];
        dudzstar_loc[ii] = (std::get<2>(du_star))[q->gt (ii)];
      }
    else
      {
        int np = q->num_parents(ii);
        for (int pp = 0; pp < np; ++pp)
          {
            dudxstar_loc[ii] += (std::get<0>(du_star))[q->gparent (pp, ii)];
            dudystar_loc[ii] += (std::get<1>(du_star))[q->gparent (pp, ii)];
            dudzstar_loc[ii] += (std::get<2>(du_star))[q->gparent (pp, ii)];
          }
        dudxstar_loc[ii] /= np;
        dudystar_loc[ii] /= np;
        dudzstar_loc[ii] /= np;
      }

  auto fun =
    [&x, &y, &z, &dudxstar_loc, &dudystar_loc, &dudzstar_loc,
     &dudx_ex, &dudy_ex, &dudz_ex]
    (double X, double Y, double Z) -> double
    {
      return
      std::pow (dudx_ex (X, Y, Z) -
                q1 (X, Y, Z, x, y, z, dudxstar_loc), 2) +
      std::pow (dudy_ex (X, Y, Z) -
                q1 (X, Y, Z, x, y, z, dudystar_loc), 2) +
      std::pow (dudz_ex (X, Y, Z) -
                q1 (X, Y, Z, x, y, z, dudzstar_loc), 2);
    };

  return std::sqrt (quad_integral (x, y, z, fun));
}


// Explicit instantiation of template functions

template
void
bim3a_advection_diffusion (tmesh_3d&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           sparse_matrix&, bool,
                           const ordering&, const ordering&);

template
void
bim3a_advection_diffusion (tmesh_3d&,
                           const std::vector<double>&,
                           const distributed_vector&,
                           sparse_matrix&, bool,
                           const ordering&, const ordering&);

/* ---- */
template
void
bim3a_reaction (tmesh_3d&,
                const std::vector<double>&,
                const std::vector<double>&,
                sparse_matrix&,
                const ordering&,
                const ordering&);

template
void
bim3a_reaction (tmesh_3d&,
                const std::vector<double>&,
                const distributed_vector&,
                sparse_matrix&,
                const ordering&,
                const ordering&);

/* ---- */
template
void
bim3a_rhs (tmesh_3d&,
           const std::vector<double>&,
           const std::vector<double>&,
           std::vector<double>&,
           const ordering&);

template
void
bim3a_rhs (tmesh_3d&,
           const std::vector<double>&,
           const distributed_vector&,
           distributed_vector&,
           const ordering&);

/* ---- */
template
void
bim3a_dirichlet_bc (tmesh_3d&, const dirichlet_bcs3&,
                    sparse_matrix&, std::vector<double>&,
                    const ordering&, const bool&);

template
void
bim3a_dirichlet_bc (tmesh_3d&, const dirichlet_bcs3&,
                    sparse_matrix&, distributed_vector&,
                    const ordering&, const bool&);

/* ---- */
template
void
bim3a_dirichlet_bc (tmesh_3d&, const dirichlet_bcs3_quad&,
                    sparse_matrix&, std::vector<double>&,
                    const ordering&, const bool&);

template
void
bim3a_dirichlet_bc (tmesh_3d&, const dirichlet_bcs3_quad&,
                    sparse_matrix&, distributed_vector&,
                    const ordering&, const bool&);

/* ---- */
template
double
nedelec_gradient (tmesh_3d::quadrant_iterator & ,
                  const std::vector<double>& , size_t);

template
double
nedelec_gradient (tmesh_3d::quadrant_iterator &,
                  const distributed_vector &, size_t);

/* ---- */
template
std::tuple<double, double, double, bool, bool, bool>
bim3c_recovered_gradient_loc (tmesh_3d::quadrant_iterator,
                              int,
                              const std::vector<double>&,
                              active_fun3);

template
std::tuple<double, double, double, bool, bool, bool>
bim3c_recovered_gradient_loc (tmesh_3d::quadrant_iterator,
                              int,
                              const distributed_vector&,
                              active_fun3);

/* ---- */
template
q2_vec3
bim3c_quadtree_pde_recovered_solution (tmesh_3d&,
                                       const std::vector<double>&,
                                       const gradient3<std::vector<double>>&);
template
q2_vec3
bim3c_quadtree_pde_recovered_solution (tmesh_3d&,
                                       const distributed_vector&,
                                       const gradient3<distributed_vector>&);
/* ---- */
template
double
estimator_grad (tmesh_3d::quadrant_iterator,
                const gradient3<std::vector<double>>&,
                const std::vector<double>&);

template
double
estimator_grad (tmesh_3d::quadrant_iterator,
                const gradient3<distributed_vector>&,
                const distributed_vector&);

/* ---- */
template
double
estimator_sol (tmesh_3d::quadrant_iterator,
               const q2_vec3 &,
               const std::vector<double>&);

template
double
estimator_sol (tmesh_3d::quadrant_iterator,
               const q2_vec3&,
               const distributed_vector&);

/* ---- */
template
double
l2_error (tmesh_3d::quadrant_iterator,
          const func3 &,
          const std::vector<double> &);

template
double
l2_error (tmesh_3d::quadrant_iterator,
          const func3 &,
          const distributed_vector &);

/* ---- */
template
double
semih1_error (tmesh_3d::quadrant_iterator,
              const func3 &,
              const func3 &,
              const func3 &,
              const std::vector<double> &);

template
double
semih1_error (tmesh_3d::quadrant_iterator,
              const func3 &,
              const func3 &,
              const func3 &,
              const distributed_vector &);

/* ---- */
template
double
semih1_star_error (tmesh_3d::quadrant_iterator,
                   const func3 &,
                   const func3 &,
                   const func3 &,
                   const gradient3<std::vector<double>> &);
template
double
semih1_star_error (tmesh_3d::quadrant_iterator,
                   const func3 &,
                   const func3 &,
                   const func3 &,
                   const gradient3<distributed_vector> &);