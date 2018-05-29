#include "quad_operators.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <set>
#include <limits>
#include <iomanip>

// Compute harmonic mean of a and b.
static double
hm (const double & a, const double & b)
{
  return 2 / (1 / a + 1 / b);
}


void
bim2a_structure (tmesh &tmsh,
                 sparse_matrix& A)
{

  std::vector<unsigned int> rows, cols;
  rows.reserve (2);
  cols.reserve (2);
  int i, j, r, c;
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
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


          for (j = 0; j < 4; ++j)
            {
              cols.clear ();
              if (! quadrant->is_hanging (j))
                cols.push_back (quadrant->gt (j));
              else
                {
                  cols.push_back (quadrant->gparent (0, j));
                  cols.push_back (quadrant->gparent (1, j));
                }
              
              for (r = 0; r < rows.size (); ++r)
                for (c = 0; c < cols.size (); ++c)
                  A[rows[r]][cols[c]] = 0.0;
            }
        }
    }
  
  A.set_properties ();
}

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
      std::array<double, 4> psi_aux;
      
      for (int n = 0; n < 4; ++n)
        {
          if (! quadrant->is_hanging (n))
            psi_aux[n] = psi[quadrant->gt (n)];
          else
            psi_aux[n] = 0.5 * (psi[quadrant->gparent (0, n)] +
                                psi[quadrant->gparent (1, n)]);
        }
      
      psi01 = psi_aux[1] - psi_aux[0];
      psi13 = psi_aux[3] - psi_aux[1];
      psi32 = psi_aux[2] - psi_aux[3];
      psi20 = psi_aux[0] - psi_aux[2];
      
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
bim2a_advection_eafe_diffusion (tmesh& mesh,
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
      std::array<double, 4> psi_aux, alpha_aux;
      
      for (int n = 0; n < 4; ++n)
        {
          if (! quadrant->is_hanging (n))
            {
              psi_aux[n] = psi[quadrant->gt (n)];
              alpha_aux[n] = alpha[quadrant->gt (n)];
            }
          else
            {
              psi_aux[n] = 0.5 * (psi[quadrant->gparent (0, n)] +
                                  psi[quadrant->gparent (1, n)]);
              alpha_aux[n] = 0.5 * (alpha[quadrant->gparent (0, n)] +
                                    alpha[quadrant->gparent (1, n)]);
            }
        }
      
      psi01 = psi_aux[1] - psi_aux[0];
      psi13 = psi_aux[3] - psi_aux[1];
      psi32 = psi_aux[2] - psi_aux[3];
      psi20 = psi_aux[0] - psi_aux[2];
      
      bimu_bernoulli(psi01, bp01, bm01);
      bimu_bernoulli(psi13, bp13, bm13);
      bimu_bernoulli(psi32, bp32, bm32);
      bimu_bernoulli(psi20, bp20, bm20);
      
      hx = quadrant->p(0, 1) - quadrant->p(0, 0);
      hy = quadrant->p(1, 2) - quadrant->p(1, 0);
      
      iel = quadrant->get_forest_quad_idx();
      
      bp01 *= hm(alpha_aux[0], alpha_aux[1]) * hy / (2 * hx);
      bm01 *= hm(alpha_aux[0], alpha_aux[1]) * hy / (2 * hx);
      bp13 *= hm(alpha_aux[1], alpha_aux[3]) * hx / (2 * hy);
      bm13 *= hm(alpha_aux[1], alpha_aux[3]) * hx / (2 * hy);
      bp32 *= hm(alpha_aux[3], alpha_aux[2]) * hy / (2 * hx);
      bm32 *= hm(alpha_aux[3], alpha_aux[2]) * hy / (2 * hx);
      bp20 *= hm(alpha_aux[2], alpha_aux[0]) * hx / (2 * hy);
      bm20 *= hm(alpha_aux[2], alpha_aux[0]) * hx / (2 * hy);
      
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
  
  double zeta_loc = 0;
  
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
            {
              rows.push_back (quadrant->gt (i));
              zeta_loc = zeta[quadrant->gt (i)];
            }
          else
            {
              rows.push_back (quadrant->gparent (0, i));
              rows.push_back (quadrant->gparent (1, i));
              zeta_loc = 0.5 * (zeta[quadrant->gparent (0, i)] +
                                zeta[quadrant->gparent (1, i)]);
            }
          
          for (int r = 0; r < rows.size (); ++r)
	    A[rows[r]][rows[r]] +=
	      (delta[iel] * zeta_loc * hx * hy / 4) /
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
   
  double g_loc = 0;
   
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
	    {
	      rows.push_back (quadrant->gt (i));
	      g_loc = g[quadrant->gt (i)];
	    }
	  else
	    {
	      rows.push_back (quadrant->gparent (0, i));
	      rows.push_back (quadrant->gparent (1, i));
	      g_loc = 0.5 * (g[quadrant->gparent (0, i)] +
			     g[quadrant->gparent (1, i)]);
	    }
            
	  for (int r = 0; r < rows.size(); ++r)
	    rhs[rows[r]] +=
	      (f[iel] * g_loc * hx * hy / 4) /
	      rows.size ();
	}
    }
}

std::vector<double>
bim2a_boundary_mass (tmesh& mesh,
		     const int & tree_idx,
		     const int & boundary_idx,
		     std::vector<double> & M,
		     const func_quad & fun)
{
  double h = 0;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
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
		  
                  M[quadrant->gt(i)] += 0.5 * h * fun (quadrant, i);
		}
	    }
	}
    }
  
  return M;
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
	      // Loop over all the boundary conditions.
              for (size_t bc = 0; bc < bcs.size (); ++bc)
                // If this boundary condition matches with
                // the current node.
                if (std::get<0> (bcs[bc]) == tree_idx
                    && std::get<1> (bcs[bc]) == boundary_idx)
                  {
		    // Mark current node so to avoid duplicate operations.
		    marked.insert (row);
		    
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
                          col = A.col_idx (j);
                          
                          if (row != col)
                            {
                              A[row][col] = 0.0;
                              
                              // If row "col" is owned by current process.
                              if (A[col].size ())
                                {
                                  rhs[col] -= A[col][row] * rhs[row];
                                  A[col][row] = 0.0;
                                }
                            }
                        }
                    
                    if (std::abs (A[row][row])
                        < std::numeric_limits<double>::epsilon ())
                      {
                        A[row][row] = std::accumulate
                          (A[row].begin (),
                           A[row].end (),
                           0.0,
                           [] (double value,
                               const std::map<int, double>::value_type & p)
                           {
                             return (value + std::abs (p.second));
                           }
                           );
                      }

		    A[row][row] *= 1e16;
                    
                    // Multiply rhs by the diagonal entry.
                    rhs[row] *= A[row][row];
                  }
            }
        }
    }
}

void
bim2a_dirichlet_bc (tmesh& mesh, const dirichlet_bcs_quad& bcs,
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
	      // Loop over all the boundary conditions.
              for (size_t bc = 0; bc < bcs.size (); ++bc)
                // If this boundary condition matches with
                // the current node.
                if (std::get<0> (bcs[bc]) == tree_idx
                    && std::get<1> (bcs[bc]) == boundary_idx)
                  {
		    // Mark current node so to avoid duplicate operations.
		    marked.insert (row); 
		    
		    // Impose boundary condition at rhs by
                    // evaluating it at the current node.
                    rhs[row] =
                      (std::get<2> (bcs[bc])) (quadrant, i);
                    
                    // Move non-diagonal entries
                    // from column "row" to rhs.
                    if (A[row].size ())
                      for (auto j = A[row].begin ();
                           j != A[row].end (); ++j)
                        {
                          col = A.col_idx (j);
                          
                          if (row != col)
                            {
                              A[row][col] = 0.0;
                              
                              // If row "col" is owned by current process.
                              if (A[col].size ())
                                {
                                  rhs[col] -= A[col][row] * rhs[row];
                                  A[col][row] = 0.0;
                                }
                            }
                        }
		    
                    if (std::abs (A[row][row])
                        < std::numeric_limits<double>::epsilon ())
                      {
                        A[row][row] = std::accumulate
                          (A[row].begin (),
                           A[row].end (),
                           0.0,
                           [] (double value,
                               const std::map<int, double>::value_type & p)
                           {
                             return (value + std::abs (p.second));
                           }
                           );
                      }

		    A[row][row] *= 1e16;
		    
                    // Multiply rhs by the diagonal entry.
                    rhs[row] *= A[row][row];
                  }
            }
        }
    }
}

// MPI_User_function.
static void replace(double *invec, double *inoutvec,
                    int *len, MPI_Datatype *dtype)
{
  for (int i = 0; i < *len; ++i)
    if (invec[i] != 0 && inoutvec[i] == 0)
      inoutvec[i] = invec[i];
}

double
nedelec_gradient (tmesh::quadrant_iterator & q,
                  const q1_vec& u, size_t i)
{
  std::array<double, 4> u_aux;
  
  double hx = q->p (0, 1) - q->p (0, 0);
  double hy = q->p (1, 2) - q->p (1, 0);
  
  for (int n = 0; n < 4; ++n)
    {
      if (! q->is_hanging (n))
        u_aux[n] = u[q->gt (n)];
      else
        u_aux[n] = 0.5 * (u[q->gparent (0, n)] +
                          u[q->gparent (1, n)]);
    }
  
  double du = 0;
  
  switch (i)
    {
    case 0:
      du = (u_aux[2] - u_aux[0]) / hy;
      break;
    case 1:
      du = (u_aux[3] - u_aux[1]) / hy;
      break;
    case 2:
      du = (u_aux[1] - u_aux[0]) / hx;
      break;
    case 3:
      du = (u_aux[3] - u_aux[2]) / hx;
      break;
    }
  
  return du;
}

gradient
bim2c_quadtree_pde_recovered_gradient (tmesh& mesh,
                                       const q1_vec& u,
                                       active_fun is_active)
{
  std::vector<double> du_x_star (mesh.num_global_nodes (), 0);
  std::vector<double> du_y_star (mesh.num_global_nodes (), 0);  
  std::vector<bool> assigned_x (mesh.num_global_nodes (), false);
  std::vector<bool> assigned_y (mesh.num_global_nodes (), false);
  
  double hx = 0, hy = 0;
  
  int node_n = 0, node_side = 0;
  std::vector<double> du_x, weights_x;
  std::vector<double> du_y, weights_y;
  
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
              (assigned_x[quadrant->gt (node)] &&
               assigned_y[quadrant->gt (node)]))
            continue;
          
          // Skip inactive quadrants.
          if (! is_active(quadrant))
            continue;
          
          du_x.clear (); weights_x.clear ();
          du_y.clear (); weights_y.clear ();
          
          du_x_star[quadrant->gt (node)] = 0;
          du_y_star[quadrant->gt (node)] = 0;
          
          // Compute Nédélec gradient on current element.
          switch (node)
            {
	    case 0:
	      du_x.push_back (nedelec_gradient(quadrant, u, 2));
	      du_y.push_back (nedelec_gradient(quadrant, u, 0));
	      break;
	    case 1:
	      du_x.push_back (nedelec_gradient(quadrant, u, 2));
	      du_y.push_back (nedelec_gradient(quadrant, u, 1));
	      break;
	    case 2:
	      du_x.push_back (nedelec_gradient(quadrant, u, 3));
	      du_y.push_back (nedelec_gradient(quadrant, u, 0));
	      break;
	    case 3:
	      du_x.push_back (nedelec_gradient(quadrant, u, 3));
	      du_y.push_back (nedelec_gradient(quadrant, u, 1));
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
			du_x.push_back (nedelec_gradient(neighbor, u, 2));
		      break;
		    case 1:
		      if (node_side == 0)
			du_x.push_back (nedelec_gradient(neighbor, u, 2));
		      break;
		    case 2:
		      if (node_side == 3)
			du_x.push_back (nedelec_gradient(neighbor, u, 3));
		      break;
		    case 3:
		      if (node_side == 2)
			du_x.push_back (nedelec_gradient(neighbor, u, 3));
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
			du_y.push_back (nedelec_gradient(neighbor, u, 0));
		      break;
		    case 1:
		      if (node_side == 3)
			du_y.push_back (nedelec_gradient(neighbor, u, 1));
		      break;
		    case 2:
		      if (node_side == 0)
			du_y.push_back (nedelec_gradient(neighbor, u, 0));
		      break;
		    case 3:
		      if (node_side == 1)
			du_y.push_back (nedelec_gradient(neighbor, u, 1));
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
              assigned_x[quadrant->gt (node)] = true;
              
              for (unsigned int ix = 0; ix < du_x.size (); ++ix)
                {
                  du_x_star[quadrant->gt (node)] +=
                    du_x[ix] * weights_x[ix];
                }
              
              du_x_star[quadrant->gt (node)] /=
                std::accumulate (weights_x.begin (),
                                 weights_x.end (), 0.0);
            }
          
          if (du_y.size () == 2)
            {
              assigned_y[quadrant->gt (node)] = true;
              
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
    }
  
  // Send data to all processes so that non-assigned values
  // on current rank get assigned by other ranks.
  std::vector<double> du_x_star_global (mesh.num_global_nodes (), 0);
  std::vector<double> du_y_star_global (mesh.num_global_nodes (), 0);
  
  MPI_Op op;
  MPI_Op_create((MPI_User_function *) replace, 1, &op);
  
  MPI_Allreduce(du_x_star.data(), du_x_star_global.data(),
                du_x_star.size(), MPI_DOUBLE,
                op, MPI_COMM_WORLD);
  
  MPI_Allreduce(du_y_star.data(), du_y_star_global.data(),
                du_y_star.size(), MPI_DOUBLE,
                op, MPI_COMM_WORLD);
  
  return std::make_pair (du_x_star_global, du_y_star_global);
}

q2_vec
bim2c_quadtree_pde_recovered_solution (tmesh& mesh,
                                       const q1_vec& u,
                                       const gradient& du)
{
  q2_vec u_star (mesh.num_local_quadrants (),
                 std::array<double, 9>({0,0,0,0,0,0,0,0,0}));
  
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
  int ix, jy, ipt;
  double wx = 0, wy = 0,
    sum = 0;
  for (ix = 0; ix < 4; ++ix)
    {
      wx = xformw (x, gw[ix]);
      
      for (jy = 0; jy < 4; ++jy)
        {
          wy = xformw (y, gw[jy]);
          sum += fun (xformx (x, gn[ix]), xformx (y, gn[jy])) * wx * wy;
        }
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
          u[4] * -8 * (X - xc) * (X - x[1]) * (Y - y[0]) * (Y - y[1]) +
          u[5] * -8 * (X - x[0]) * (X - xc) * (Y - y[0]) * (Y - y[1]) +
          u[6] * -8 * (X - x[0]) * (X - x[1]) * (Y - yc) * (Y - y[1]) +
          u[7] * -8 * (X - x[0]) * (X - x[1]) * (Y - y[0]) * (Y - yc) +
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
      if (! q->is_hanging (ii))
        {
          dudxstar_loc[ii] = (du_star.first)[q->gt(ii)];
          dudystar_loc[ii] = (du_star.second)[q->gt(ii)];
          u_loc[ii] = u[q->gt(ii)];
        }
      else
        {
          dudxstar_loc[ii] = 0.5 * ((du_star.first)[q->gparent(0, ii)] +
                                    (du_star.first)[q->gparent(1, ii)]);
          dudystar_loc[ii] = 0.5 * ((du_star.second)[q->gparent(0, ii)] +
                                    (du_star.second)[q->gparent(1, ii)]);
          u_loc[ii] = 0.5 * (u[q->gparent(0, ii)] +
                             u[q->gparent(1, ii)]);
        }
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
      if (! q->is_hanging (ii))
        u_loc[ii] = u[q->gt(ii)];
      else
        u_loc[ii] = 0.5 * (u[q->gparent(0, ii)] +
                           u[q->gparent(1, ii)]);
    }

  auto fun =
    [x, y, ustar_loc, u_loc]
    (double X, double Y) -> double
    {
      return
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

// Compute ||u - u_ex||_L^2(q).
double
l2_error (tmesh::quadrant_iterator q,
          const func & u_ex,
          const q1_vec & u)
{
  double
    x[2] = {q->p(0,0), q->p(0,1)},
    y[2] = {q->p(1,0), q->p(1,3)};

  double u_loc[4] = {0,0,0,0};

  for (int ii = 0; ii < 4; ++ii)
    {
      if (! q->is_hanging (ii))
        u_loc[ii] = u[q->gt(ii)];
      else
        u_loc[ii] = 0.5 * (u[q->gparent(0, ii)] +
                           u[q->gparent(1, ii)]);
    }

  auto fun =
    [x, y, u_loc, u_ex]
    (double X, double Y) -> double
    {
      return
      std::pow (q1 (X, Y, x, y, u_loc) - u_ex(X, Y), 2);
    };
    
  return std::sqrt(quad_integral (x, y, fun));
}

// Compute |u - u_ex|_H^1(q).
double
semih1_error (tmesh::quadrant_iterator q,
              const func & dudx_ex,
              const func & dudy_ex,
              const q1_vec & u)
{
  double
    x[2] = {q->p(0,0), q->p(0,1)},
    y[2] = {q->p(1,0), q->p(1,3)};

  double u_loc[4] = {0,0,0,0};

  for (int ii = 0; ii < 4; ++ii)
    {
      if (! q->is_hanging (ii))
        u_loc[ii] = u[q->gt(ii)];
      else
        u_loc[ii] = 0.5 * (u[q->gparent(0, ii)] +
                           u[q->gparent(1, ii)]);
    }

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
    
  return std::sqrt(quad_integral (x, y, fun));
}

// Compute ||u_star - u_ex||_L^2(q).
double
l2_star_error (tmesh::quadrant_iterator q,
               const func & u_ex,
               const q2_vec & ustar)
{
  double
    x[2] = {q->p(0,0), q->p(0,1)},
    y[2] = {q->p(1,0), q->p(1,3)};
  
  double ustar_loc[9] = {0,0,0,0,0,0,0,0,0};
  
  for (int ii = 0; ii < 9; ++ii)
    {
      ustar_loc[ii] = (ustar[q->get_forest_quad_idx()])[ii];
    }
  
  auto fun =
    [x, y, ustar_loc, u_ex]
    (double X, double Y) -> double
    {
      return
      std::pow (q2 (X, Y, x, y, ustar_loc) - u_ex(X, Y), 2);
    };
  
  return std::sqrt(quad_integral (x, y, fun));
}


// Compute ||du_star - grad(u_ex)||_L^2(q).
double
semih1_star_error (tmesh::quadrant_iterator q,
                   const func & dudx_ex,
                   const func & dudy_ex,
                   const gradient & du_star)
{
  double
    x[2] = {q->p(0,0), q->p(0,1)},
    y[2] = {q->p(1,0), q->p(1,3)};
  
  double dudxstar_loc[4] = {0,0,0,0};
  double dudystar_loc[4] = {0,0,0,0};
  
  for (int ii = 0; ii < 4; ++ii)
    {
      if (! q->is_hanging (ii))
        {
          dudxstar_loc[ii] = (du_star.first)[q->gt(ii)];
          dudystar_loc[ii] = (du_star.second)[q->gt(ii)];
        }
      else
        {
          dudxstar_loc[ii] = 0.5 * ((du_star.first)[q->gparent(0, ii)] +
                                    (du_star.first)[q->gparent(1, ii)]);
          dudystar_loc[ii] = 0.5 * ((du_star.second)[q->gparent(0, ii)] +
                                    (du_star.second)[q->gparent(1, ii)]);
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
  
  return std::sqrt(quad_integral (x, y, fun));
}
