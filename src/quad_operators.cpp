#include "quad_operators.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <set>

void 
bim2a_advection_diffusion (tmesh& mesh,
                           const std::vector<double>& alpha,
                           const std::vector<double>& psi,
                           sparse_matrix& A)
{
  double psi01 = 0;
  double psi13 = 0;
  double psi32 = 0;
  double psi20 = 0;
  
  double bp01 = 0, bm01 = 0;
  double bp13 = 0, bm13 = 0;
  double bp32 = 0, bm32 = 0;
  double bp20 = 0, bm20 = 0;
  
  double hx = 0, hy = 0;
  
  std::array<std::array<double, 4>, 4> Aloc;
  
  unsigned int iel = 0;
  std::vector<unsigned int> rows, cols;
  rows.reserve(2);
  cols.reserve(2);
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      psi01 = psi[quadrant->t(1)] - psi[quadrant->t(0)];
      psi13 = psi[quadrant->t(3)] - psi[quadrant->t(1)];
      psi32 = psi[quadrant->t(2)] - psi[quadrant->t(3)];
      psi20 = psi[quadrant->t(0)] - psi[quadrant->t(2)];
      
      bimu_bernoulli(psi01, bp01, bm01);
      bimu_bernoulli(psi13, bp13, bm13);
      bimu_bernoulli(psi32, bp32, bm32);
      bimu_bernoulli(psi20, bp20, bm20);
      
      hx = quadrant->p(0, 1) - quadrant->p(0, 0);
      hy = quadrant->p(1, 2) - quadrant->p(1, 0);
      
      iel = quadrant->get_forest_quad_idx();
      
      bp01 *= alpha[iel] * hy / (2 * hx);
      bm01 *= alpha[iel] * hy / (2 * hx);
      bp13 *= alpha[iel] * hx / (2 * hy);
      bm13 *= alpha[iel] * hx / (2 * hy);
      bp32 *= alpha[iel] * hy / (2 * hx);
      bm32 *= alpha[iel] * hy / (2 * hx);
      bp20 *= alpha[iel] * hx / (2 * hy);
      bm20 *= alpha[iel] * hx / (2 * hy);
      
      Aloc[0] = { bm01 + bp20, -bp01,        -bm20,         0          };
      Aloc[1] = {-bm01,         bp01 + bm13,  0,           -bp13       };
      Aloc[2] = {-bp20,         0,            bp32 + bm20, -bm32       };
      Aloc[3] = { 0,           -bm13,        -bp32,         bm32 + bp13};
      
      for(int i = 0; i < 4; ++i)
        {
          rows.clear();
          
          if (!quadrant->is_hanging(i))
            rows.push_back (quadrant->gt(i));
          else
            {
              rows.push_back (quadrant->gparent(0, i));
              rows.push_back (quadrant->gparent(1, i));
            }
          
          for(int j = 0; j < 4; ++j)
            {
              cols.clear();
              
              if (!quadrant->is_hanging(j))
                cols.push_back (quadrant->gt(j));
              else
                {
                  cols.push_back (quadrant->gparent(0, j));
                  cols.push_back (quadrant->gparent(1, j));
                }
              
              for (int r = 0; r < rows.size(); ++r)
                for (int c = 0; c < cols.size(); ++c)
                  {
                    A[rows[r]][cols[c]] += Aloc[i][j] /
                      (rows.size() * cols.size());
                  }
            }
        }
    }
}

void
bim2a_reaction (tmesh& mesh,
                const std::vector<double>& delta,
                const std::vector<double>& zeta,
                sparse_matrix& A)
{
  double hx = 0, hy = 0;
  
  unsigned int iel = 0;
  std::vector<unsigned int> rows;
  rows.reserve(2);
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p(0, 1) - quadrant->p (0, 0);
      hy = quadrant->p(1, 2) - quadrant->p (1, 0);
      
      iel = quadrant->get_forest_quad_idx ();
      
      for(int i = 0; i < 4; ++i)
        {
          rows.clear();
          
          if (!quadrant->is_hanging (i))
            rows.push_back (quadrant->gt (i));
          else
            {
              rows.push_back (quadrant->gparent (0, i));
              rows.push_back (quadrant->gparent (1, i));
            }
          
          for (int r = 0; r < rows.size (); ++r)
            A[rows[r]][rows[r]] +=
              (delta[iel] * zeta[quadrant->t (i)] * hx * hy / 4) /
              rows.size ();
        }
    }
}

void
bim2a_rhs (tmesh& mesh,
           const std::vector<double>& f,
           const std::vector<double>& g,
           std::vector<double>& rhs)
{
   double hx = 0, hy = 0;
   
   unsigned int iel = 0;
   std::vector<unsigned int> rows;
   rows.reserve (2);
   
   for (auto quadrant = mesh.begin_quadrant_sweep ();
        quadrant != mesh.end_quadrant_sweep ();
        ++quadrant)
     {
        hx = quadrant->p (0, 1) - quadrant->p (0, 0);
        hy = quadrant->p (1, 2) - quadrant->p (1, 0);
        
        iel = quadrant->get_forest_quad_idx ();
        
        for(int i = 0; i < 4; ++i)
          {
            rows.clear ();
            
            if (! quadrant->is_hanging (i))
              rows.push_back (quadrant->gt (i));
            else
              {
                rows.push_back (quadrant->gparent (0, i));
                rows.push_back (quadrant->gparent (1, i));
              }
            
            for (int r = 0; r < rows.size(); ++r)
              rhs[rows[r]] +=
                (f[iel] * g[quadrant->t (i)] * hx * hy / 4) /
                rows.size ();
          }
     }
}

void
bim2a_dirichlet_bc (tmesh& mesh, const dirichlet_bcs& bcs,
                    sparse_matrix& A, std::vector<double>& rhs)
{
  int boundary_idx, tree_idx;
  unsigned int row, col;
  
  std::set<unsigned int> marked;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      tree_idx = quadrant->get_tree_idx ();
      
      for (int i = 0; i < 4; ++i)
        {
          boundary_idx = quadrant->e (i);
          row = quadrant->gt (i);
          
          // If current node is on boundary and has not
          // been handled before.
          if (boundary_idx != tmesh::quadrant_t::NOT_ON_BOUNDARY
              && marked.count(row) == 0)
            {
              // Mark current node so to avoid duplicate operations.
              marked.insert (row); 
              
              // Loop over all the boundary conditions.
              for (size_t bc = 0; bc < bcs.size (); ++bc)
                // If this boundary condition matches with
                // the current node.
                if (std::get<0> (bcs[bc]) == tree_idx
                    && std::get<1> (bcs[bc]) == boundary_idx)
                  {
                    // Impose boundary condition at rhs by
                    // evaluating it at the current node.
                    rhs[row] =
                      (std::get<2> (bcs[bc]))
                      (quadrant->p (0, i), quadrant->p (1, i));
                      
                    // Move non-diagonal entries
                    // from column "row" to rhs.
                    if (A[row].size ())
                      for (auto j = A[row].begin ();
                           j != A[row].end (); ++j)
                        {
                          col = A.col_idx(j);
                              
                          if (row != col)
                            {
                              A[row][col] = 0.0;
                              rhs[col] -= A[col][row] * rhs[row];
                              A[col][row] = 0.0;
                            }
                        }
                        
                      
                    // Multiply rhs by the diagonal entry.
                    rhs[row] *= A[row][row];
                  }
            }
        }
    }
}

gradient
bim2c_quadtree_pde_recovered_gradient (tmesh& mesh, const q1_vec& u)
{
  std::vector<double> du_x_star (mesh.num_global_nodes (), 0);
  std::vector<double> du_y_star (mesh.num_global_nodes (), 0);  
  std::vector<bool> assigned (mesh.num_global_nodes (), false);
  
  double hx = 0, hy = 0;
  
  int node_n = 0, node_side = 0;
  std::vector<double> du_x, weights_x;
  std::vector<double> du_y, weights_y;
  std::vector<double> u_aux;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      // Loop over vertices of current quadrant that
      // are non-hanging and have not been processed yet.
      for (int node = 0; node < 4; ++node)
        {
          hx = quadrant->p (0, 1) - quadrant->p (0, 0);
          hy = quadrant->p (1, 2) - quadrant->p (1, 0);
      
          if (quadrant->is_hanging (node) ||
              assigned[quadrant->gt (node)])
            continue;
          
          du_x.clear (); weights_x.clear ();
          du_y.clear (); weights_y.clear ();
          u_aux.clear (); u_aux.resize (4);
          
          // Compute Nédélec gradient on current element.
          for (int n = 0; n < 4; ++n)
            {
              if (! quadrant->is_hanging (n))
                u_aux[n] = u[quadrant->gt (n)];
              else
                u_aux[n] = 0.5 * (u[quadrant->gparent (0, n)] +
                                  u[quadrant->gparent (1, n)]);
            }
          
          switch (node)
            {
              case 0:
                du_x.push_back ((u_aux[1] - u[quadrant->gt (0)]) / hx);
                du_y.push_back ((u_aux[2] - u[quadrant->gt (0)]) / hy);
                break;
              case 1:
                du_x.push_back ((u[quadrant->gt (1)] - u_aux[0]) / hx);
                du_y.push_back ((u_aux[3] - u[quadrant->gt (1)]) / hy);
                break;
              case 2:
                du_x.push_back ((u_aux[3] - u[quadrant->gt (2)]) / hx);
                du_y.push_back ((u[quadrant->gt (2)] - u_aux[0]) / hy);
                break;
              case 3:
                du_x.push_back((u[quadrant->gt (3)] - u_aux[2]) / hx);
                du_y.push_back((u[quadrant->gt (3)] - u_aux[1]) / hy);
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
                      du_x.push_back ((u[neighbor->gt (1)] -
                                       u[neighbor->gt (0)]) / hx);
                    if (node == 2)
                      du_y.push_back ((u[neighbor->gt (2)] -
                                       u[neighbor->gt (0)]) / hy);
                    break;
                  case 1:
                    if (node == 0)
                      du_x.push_back ((u[neighbor->gt (1)] -
                                       u[neighbor->gt (0)]) / hx);
                    if (node == 3)
                      du_y.push_back ((u[neighbor->gt (3)] -
                                       u[neighbor->gt (1)]) / hy);
                    break;
                  case 2:
                    if (node == 3)
                      du_x.push_back ((u[neighbor->gt (3)] -
                                       u[neighbor->gt (2)]) / hx);
                    if (node == 0)
                      du_y.push_back ((u[neighbor->gt (2)] -
                                       u[neighbor->gt (0)]) / hy);
                    break;
                  case 3:
                    if (node == 2)
                      du_x.push_back ((u[neighbor->gt (3)] -
                                       u[neighbor->gt (2)]) / hx);
                    if (node == 1)
                      du_y.push_back ((u[neighbor->gt (3)] -
                                       u[neighbor->gt (1)]) / hy);
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
                          du_x.push_back ((u[neighbor->gt (1)] -
                                           u[neighbor->gt (0)]) / hx);
                        break;
                      case 1:
                        if (node_side == 0)
                          du_x.push_back ((u[neighbor->gt (1)] -
                                           u[neighbor->gt (0)]) / hx);
                        break;
                      case 2:
                        if (node_side == 3)
                          du_x.push_back ((u[neighbor->gt (3)] -
                                           u[neighbor->gt (2)]) / hx);
                        break;
                      case 3:
                        if (node_side == 2)
                          du_x.push_back ((u[neighbor->gt (3)] -
                                           u[neighbor->gt (2)]) / hx);
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
                          du_y.push_back((u[neighbor->gt (2)] -
                                          u[neighbor->gt (0)]) / hy);
                        break;
                      case 1:
                        if (node_side == 3)
                          du_y.push_back ((u[neighbor->gt (3)] -
                                           u[neighbor->gt (1)]) / hy);
                        break;
                      case 2:
                        if (node_side == 0)
                          du_y.push_back ((u[neighbor->gt (2)] -
                                           u[neighbor->gt (0)]) / hy);
                        break;
                      case 3:
                        if (node_side == 1)
                          du_y.push_back ((u[neighbor->gt (3)] -
                                           u[neighbor->gt (1)]) / hy);
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
          
          if (du_x.size () < 2 || du_y.size () < 2)
            continue;
          else
            assigned[quadrant->gt (node)] = true;
          
          for (unsigned int ix = 0; ix < du_x.size (); ++ix)
            {
              du_x_star[quadrant->gt (node)] +=
                du_x[ix] * weights_x[ix];
            }
          
          du_x_star[quadrant->gt (node)] /=
            std::accumulate (weights_x.begin (),
                             weights_x.end (), 0.0);
          
          for (unsigned int iy = 0; iy < du_y.size (); ++iy)
            {
              du_y_star[quadrant->gt (node)] +=
                du_y[iy] * weights_y[iy];
            }
          
          du_y_star[quadrant->gt (node)] /=
            std::accumulate (weights_y.begin (),
                             weights_y.end (), 0.0);
        }
    }
  
  return std::make_pair (du_x_star, du_y_star);
}

q2_vec
bim2c_quadtree_pde_recovered_solution (tmesh& mesh,
                                       const q1_vec& u,
                                       const gradient& du)
{
  q2_vec u_star (mesh.num_local_quadrants (),
                 std::array<double, 9>({0,0,0,0,0,0,0,0,0}));
  
  double hx = 0, hy = 0;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p (0, 1) - quadrant->p (0, 0);
      hy = quadrant->p (1, 2) - quadrant->p (1, 0);
      
      // Compute values at vertices.
      for (int i = 0; i < 4; ++i)
        u_star[quadrant->get_forest_quad_idx ()][i] =
          u[quadrant->gt(i)];
      
      // Compute values at faces.
      u_star[quadrant->get_forest_quad_idx ()][4] =
        0.5 * (u[quadrant->gt(0)] + u[quadrant->gt(1)])
        + hx * (du.first[quadrant->gt(0)]
                - du.first[quadrant->gt(1)]) / 8;
      
      u_star[quadrant->get_forest_quad_idx ()][5] =
        0.5 * (u[quadrant->gt(2)] + u[quadrant->gt(3)])
        + hx * (du.first[quadrant->gt(2)]
                - du.first[quadrant->gt(3)]) / 8;
      
      u_star[quadrant->get_forest_quad_idx ()][6] =
        0.5 * (u[quadrant->gt(0)] + u[quadrant->gt(2)])
        + hy * (du.second[quadrant->gt(0)]
                - du.second[quadrant->gt(2)]) / 8;
      
      u_star[quadrant->get_forest_quad_idx ()][7] =
        0.5 * (u[quadrant->gt(1)] + u[quadrant->gt(3)])
        + hy * (du.second[quadrant->gt(1)]
                - du.second[quadrant->gt(3)]) / 8;
      
      // Compute value at cell midpoint.
      u_star[quadrant->get_forest_quad_idx ()][8] =
        0.25 * (u_star[quadrant->get_forest_quad_idx ()][4] +
                u_star[quadrant->get_forest_quad_idx ()][5] +
                u_star[quadrant->get_forest_quad_idx ()][6] +
                u_star[quadrant->get_forest_quad_idx ()][7])
        + hx * 0.5 * (du.first[quadrant->gt(0)]
                      + du.first[quadrant->gt(2)]
                      - du.first[quadrant->gt(1)]
                      - du.first[quadrant->gt(3)]) / 16
        + hy * 0.5 * (du.second[quadrant->gt(0)]
                      + du.second[quadrant->gt(1)]
                      - du.second[quadrant->gt(2)]
                      - du.second[quadrant->gt(3)]) / 16;
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
  int ix, jy, ipt;
  double wx = 0,
    sum = 0;
  for (ix = 0; ix < 4; ++ix)
    {
      wx = xformw (x, gw[ix]);
      for (jy = 0; jy < 4; ++jy)
        sum += fun (xformx (x, gn[ix]), xformx (y, gn[jy])) *
          wx * xformw (y, gw[jy]);
    }
  return (sum);
}

// Evaluate Nédelec x-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]])
// at (X, Y).
static double
dudx (double X, double Y, const double *x,
      const double *y, const double *u)
{
  double hx = (x[1] - x[0]);
  double hy = (y[1] - y[0]);
  
  double db = (u[1] - u[0]) / hx;
  double dt = (u[3] - u[2]) / hx;
  
  return (db * (y[1] - Y) + dt * (Y - y[0])) / hy;
}

// Evaluate Nédelec y-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]])
// at (X, Y).
static double
dudy (double X, double Y, const double *x,
      const double *y, const double *u)
{
  double hx = (x[1] - x[0]);
  double hy = (y[1] - y[0]);
  
  double dl = (u[2] - u[0]) / hy;
  double dr = (u[3] - u[1]) / hy;
  
  return (dl * (x[1] - X) + dr * (X - x[0])) / hx;
}

// Evaluate u (using Q1 basis functions
// on quadrant [x[0], x[1]] x [y[0], y[1]])
// at (X, Y).
static double
q1 (double X, double Y, const double *x,
    const double *y, const double *u)
{
  double hx = (x[1] - x[0]);
  double hy = (y[1] - y[0]);

  return ((u[0] * (X - x[1]) * (Y - y[1]) +
           u[1] * -(X - x[0]) * (Y - y[1]) +
           u[2] * -(X - x[1]) * (Y - y[0]) +
           u[3] * (X - x[0]) * (Y - y[0])) /
          (hx * hy));
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

  double hx = (x[1] - x[0]);
  double hy = (y[1] - y[0]);

  return (u[0] * 4 * (X - xc) * (X - x[1]) * (Y - yc) * (Y - y[1]) +
          u[1] * 4 * (X - x[0]) * (X - xc) * (Y - yc) * (Y - y[1]) +
          u[2] * 4 * (X - xc) * (X - x[1]) * (Y - y[0]) * (Y - yc) +
          u[3] * 4 * (X - x[0]) * (X - xc) * (Y - y[0]) * (Y - yc) +
          u[4] * -8 * (X - x[0]) * (X - x[1]) * (Y - yc) * (Y - y[1]) +
          u[5] * -8 * (X - x[0]) * (X - x[1]) * (Y - y[0]) * (Y - yc) +
          u[6] * -8 * (X - xc) * (X - x[1]) * (Y - y[0]) * (Y - y[1]) +
          u[7] * -8 * (X - x[0]) * (X - xc) * (Y - y[0]) * (Y - y[1]) +
          u[8] * 16 * (X - x[0]) * (X - x[1]) * (Y - y[0]) * (Y - y[1])) /
         (hx * hx * hy * hy);
}

// Compute ||grad^* u - grad u||_L^2(q).
double estimator_grad(tmesh::quadrant_iterator q,
                      const gradient & du_star,
                      const q1_vec & u)
{
  double
    x[2] = {q->p(0,0), q->p(0,1)},
    y[2] = {q->p(1,0), q->p(1,3)};
  
  double dudxstar_loc[4] = {0,0,0,0};
  double dudystar_loc[4] = {0,0,0,0};
  double u_loc[4] = {0,0,0,0};
  
  for (int ii = 0; ii < 4; ++ii)
    {
      dudxstar_loc[ii] = (du_star.first)[q->gt(ii)];
      dudystar_loc[ii] = (du_star.second)[q->gt(ii)];
      u_loc[ii] = u[q->gt(ii)];
    }
  
  auto fun =
    [x, y, dudxstar_loc, dudystar_loc, u_loc]
    (double X, double Y) -> double
    {
      double err =
      std::pow (dudy (X, Y, x, y, u_loc) -
                q1 (X, Y, x, y, dudystar_loc), 2) +
      std::pow (dudx (X, Y, x, y, u_loc) -
                q1 (X, Y, x, y, dudxstar_loc), 2);
    };
  
  return std::sqrt(quad_integral (x, y, fun));
}

// Refinement marker function based on ZZ estimator
// for the recovered gradient du*.
int 
zz_marker_grad (tmesh::quadrant_iterator q,
                const gradient & du_star,
                const q1_vec & u,
                double limit)
{
  return estimator_grad(q, du_star, u) > limit ? 1 : 0;
}

// Compute ||u^* - u||_L^2(q).
double estimator_sol(tmesh::quadrant_iterator q,
                     const q2_vec & ustar,
                     const q1_vec & u)
{
  double
    x[2] = {q->p(0,0), q->p(0,1)},
    y[2] = {q->p(1,0), q->p(1,3)};

  double ustar_loc[9] = {0,0,0,0,0,0,0,0,0};
  double u_loc[4] = {0,0,0,0};

  for (int ii = 0; ii < 9; ++ii)
    {
      ustar_loc[ii] = (ustar[q->get_forest_quad_idx()])[ii];
    }
  
  for (int ii = 0; ii < 4; ++ii)
    {
      u_loc[ii] = u[q->gt(ii)];
    }

  auto fun =
    [x, y, ustar_loc, u_loc]
    (double X, double Y) -> double
    {
      double err =
      std::pow (q1 (X, Y, x, y, u_loc) -
                q2 (X, Y, x, y, ustar_loc), 2);
    };
    
  return std::sqrt(quad_integral (x, y, fun));
}

// Refinement marker function based on ZZ estimator
// for the recovered solution u*.
int 
zz_marker_sol (tmesh::quadrant_iterator q,
               const q2_vec & ustar,
               const q1_vec & u,
               double limit)
{
  return estimator_sol(q, ustar, u) > limit ? 1 : 0;
}

