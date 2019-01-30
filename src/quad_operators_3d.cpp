#include "quad_operators_3d.h"

#include <cmath>
#include <functional>
#include <numeric>
#include <set>
#include <limits>
#include <iomanip>

// MPI_User_function.
static void replace(double *invec, double *inoutvec,
                    int *len, MPI_Datatype *dtype)
{
  for (int i = 0; i < *len; ++i)
    if (invec[i] != 0 && inoutvec[i] == 0)
      inoutvec[i] = invec[i];
}

void
bim3a_structure (tmesh_3d &tmsh,
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
      for(int i = 0; i < 8; ++i)
        {
          rows.clear();
          
          if (!quadrant->is_hanging(i))
            rows.push_back (quadrant->gt(i));
          else
            for (int pp = 0; pp < quadrant->num_parents (i); ++pp)
              rows.push_back (quadrant->gparent(pp, i));
          
          for(int j = 0; j < 8; ++j)
            {
              cols.clear();
              
              if (!quadrant->is_hanging(j))
                cols.push_back (quadrant->gt(j));
              else
                for (int pp = 0; pp < quadrant->num_parents (j); ++pp)
                  cols.push_back (quadrant->gparent(pp, j));
              
              for (int r = 0; r < rows.size(); ++r)
                for (int c = 0; c < cols.size(); ++c)
                  {
                    A[rows[r]][cols[c]] = 0.0;
                  }
            }
        }
    }
  
  A.set_properties ();
}

void
bim3a_advection_diffusion (tmesh_3d& mesh,
                           const std::vector<double>& alpha,
                           const std::vector<double>& psi,
                           sparse_matrix& A)
{
  
  double bp01, bp13, bp23, bp02, bp04, bp15,
    bp26, bp37, bp45, bp57, bp67 ,bp46;
  
  double bm01, bm13, bm23, bm02, bm04, bm15,
    bm26, bm37, bm45, bm57, bm67 ,bm46;
  
  double hx, hy, hz;
  
  std::array<std::array<double, 8>, 8> Aloc;
  
  unsigned int iel = 0;
  std::vector<unsigned int> rows, cols;
  rows.reserve(4);
  cols.reserve(4);
  
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
      
      for(int i = 0; i < 8; ++i)
        {
          rows.clear();
          
          if (!quadrant->is_hanging(i))
            rows.push_back (quadrant->gt(i));
          else
            for (int pp = 0; pp < quadrant->num_parents (i); ++pp)
              rows.push_back (quadrant->gparent(pp, i));
          
          for(int j = 0; j < 8; ++j)
            {
              cols.clear();
              
              if (!quadrant->is_hanging(j))
                cols.push_back (quadrant->gt(j));
              else
                for (int pp = 0; pp < quadrant->num_parents (j); ++pp)
                  cols.push_back (quadrant->gparent(pp, j));
              
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
bim3a_reaction (tmesh_3d& mesh,
                const std::vector<double>& delta,
                const std::vector<double>& zeta,
                sparse_matrix& A)
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
            A[rows[r]][rows[r]] +=
              (delta[iel] * z_loc * hx * hy * hz / 8) / rows.size ();
        }
    }
}

void
bim3a_rhs (tmesh_3d& mesh,
           const std::vector<double>& f,
           const std::vector<double>& g,
           std::vector<double>& rhs)
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

std::vector<double>
bim3a_boundary_mass (tmesh_3d & mesh,
		     const int & tree_idx,
		     const int & boundary_idx,
		     std::vector<double> & M,
		     const func3_quad & fun)
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

void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3& bcs,
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
		    
		    // Impose boundary condition at rhs by
                    // evaluating it at the current node.
                    rhs[row] =
                      (std::get<2> (bcs[bc]))
                      (quadrant->p (0, i),
                       quadrant->p (1, i),
                       quadrant->p (2, i));
                    
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
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs3_quad& bcs,
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

/* CCI: BEGIN ADDED */

/// Edge ordering:
///
///             ______7__________     
///            /|                /|
///           4 |               5 |
///          /  |              /  |
///         /___|____6________/   |
///        |    |             |   11
///        |    10            |   |
///        |    |             9   |
///        8    |_______3_____|___|
///        |   /              |   /
///        |  0               |  1
///        | /                | /
///        |/_______2_________|/
///
///

static constexpr
std::array<std::array<int, 3>, 12> edge =
  {0,2,1, 1,3,1, 0,1,0, 2,3,0, 4,6,1, 5,7,1, 4,5,0, 7,6,0,
    0,4,2, 1,5,2, 2,6,2, 3,7,2};

double
nedelec_gradient (tmesh_3d::quadrant_iterator & q,
                  const q1_vec& u, size_t i)
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

std::tuple<double, double, double, bool, bool, bool>
bim2c_recovered_gradient_loc (tmesh_3d::quadrant_iterator quadrant,
                              int node,
                              const q1_vec& u,
                              active_fun3 is_active)
{
  double hx = quadrant->p (0, 1) - quadrant->p (0, 0);
  double hy = quadrant->p (1, 2) - quadrant->p (1, 0);
  double hz = quadrant->p (2, 4) - quadrant->p (2, 0);

  int node_n = 0, node_side = 0;
  std::vector<double> du_x, weights_x;
  std::vector<double> du_y, weights_y;
  std::vector<double> du_z, weights_z;

  du_x.clear (); weights_x.clear ();
  du_y.clear (); weights_y.clear ();
  du_z.clear (); weights_z.clear ();

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
      du_x.push_back (nedelec_gradient (quadrant, u, 2));
      du_y.push_back (nedelec_gradient (quadrant, u, 0));
      du_z.push_back (nedelec_gradient (quadrant, u, 8));
      break;
    case 1:
      du_x.push_back (nedelec_gradient (quadrant, u, 2));
      du_y.push_back (nedelec_gradient (quadrant, u, 1));
      du_z.push_back (nedelec_gradient (quadrant, u, 9));
      break;
    case 2:
      du_x.push_back (nedelec_gradient (quadrant, u, 3));
      du_y.push_back (nedelec_gradient (quadrant, u, 0));
      du_z.push_back (nedelec_gradient (quadrant, u, 10));
      break;
    case 3:
      du_x.push_back (nedelec_gradient (quadrant, u, 3));
      du_y.push_back (nedelec_gradient (quadrant, u, 1));
      du_z.push_back (nedelec_gradient (quadrant, u, 11));
      break;
    case 4:
      du_x.push_back (nedelec_gradient (quadrant, u, 6));
      du_y.push_back (nedelec_gradient (quadrant, u, 4));
      du_z.push_back (nedelec_gradient (quadrant, u, 8));
      break;
    case 5:
      du_x.push_back (nedelec_gradient (quadrant, u, 6));
      du_y.push_back (nedelec_gradient (quadrant, u, 5));
      du_z.push_back (nedelec_gradient (quadrant, u, 9));
      break;
    case 6:
      du_x.push_back (nedelec_gradient (quadrant, u, 7));
      du_y.push_back (nedelec_gradient (quadrant, u, 4));
      du_z.push_back (nedelec_gradient (quadrant, u, 10));
      break;
    case 7:
      du_x.push_back (nedelec_gradient (quadrant, u, 7));
      du_y.push_back (nedelec_gradient (quadrant, u, 5));
      du_z.push_back (nedelec_gradient (quadrant, u, 11));
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

      hx = neighbor->p (0, 1) - neighbor->p (0, 0);
      hy = neighbor->p (1, 2) - neighbor->p (1, 0);
      hz = neighbor->p (2, 4) - neighbor->p (2, 0);

      switch (node_n)
        {
        case 0:
          if (node == 1)
            du_x.push_back (nedelec_gradient(neighbor, u, 2));
          if (node == 2)
            du_y.push_back (nedelec_gradient(neighbor, u, 0));
          if (node == 4)
            du_z.push_back (nedelec_gradient(neighbor, u, 8));
          break;
        case 1:
          if (node == 0)
            du_x.push_back (nedelec_gradient(neighbor, u, 2));
          if (node == 3)
            du_y.push_back (nedelec_gradient(neighbor, u, 1));
          if (node == 5)
            du_z.push_back (nedelec_gradient(neighbor, u, 9));
          break;
        case 2:
          if (node == 3)
            du_x.push_back (nedelec_gradient(neighbor, u, 3));
          if (node == 0)
            du_y.push_back (nedelec_gradient(neighbor, u, 0));
          if (node == 6)
            du_z.push_back (nedelec_gradient(neighbor, u, 10));
          break;
        case 3:
          if (node == 2)
            du_x.push_back (nedelec_gradient(neighbor, u, 3));
          if (node == 1)
            du_y.push_back (nedelec_gradient(neighbor, u, 1));
          if (node == 7)
            du_z.push_back (nedelec_gradient(neighbor, u, 11));
          break;
        case 4:
          if (node == 5)
            du_x.push_back (nedelec_gradient(neighbor, u, 6));
          if (node == 6)
            du_y.push_back (nedelec_gradient(neighbor, u, 4));
          if (node == 0)
            du_z.push_back (nedelec_gradient(neighbor, u, 8));
          break;
        case 5:
          if (node == 4)
            du_x.push_back (nedelec_gradient(neighbor, u, 6));
          if (node == 7)
            du_y.push_back (nedelec_gradient(neighbor, u, 5));
          if (node == 1)
            du_z.push_back (nedelec_gradient(neighbor, u, 9));
          break;
        case 6:
          if (node == 7)
            du_x.push_back (nedelec_gradient(neighbor, u, 7));
          if (node == 4)
            du_y.push_back (nedelec_gradient(neighbor, u, 4));
          if (node == 2)
            du_z.push_back (nedelec_gradient(neighbor, u, 10));
          break;
        case 7:
          if (node == 6)
            du_x.push_back (nedelec_gradient(neighbor, u, 7));
          if (node == 5)
            du_y.push_back (nedelec_gradient(neighbor, u, 5));
          if (node == 3)
            du_z.push_back (nedelec_gradient(neighbor, u, 11));
          break;
        }

      if (weights_x.size () < du_x.size ())
        weights_x.push_back (1 / hx);

      if (weights_y.size () < du_y.size ())
        weights_y.push_back (1 / hy);

      if (weights_z.size () < du_z.size ())
        weights_z.push_back (1 / hz);    
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

          for (node_n = 0; node_n < 8; ++node_n)
            if (neighbor->gt (node_n) ==
                quadrant->gt (node_side))
              break;

          // If not, or if node_n is hanging,
          // switch to the next neighbor.
          if (node_n == 8 || neighbor->is_hanging (node_n))
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
            case 4:
              if (node_side == 5)
                du_x.push_back (nedelec_gradient (neighbor, u, 6));
              break;
            case 5:
              if (node_side == 4)
                du_x.push_back (nedelec_gradient (neighbor, u, 6));
              break;
            case 6:
              if (node_side == 7)
                du_x.push_back (nedelec_gradient (neighbor, u, 7));
              break;
            case 7:
              if (node_side == 6)
                du_x.push_back (nedelec_gradient (neighbor, u, 7));
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

          for (node_n = 0; node_n < 8; ++node_n)
            if (neighbor->gt (node_n) ==
                quadrant->gt (node_side))
              break;

          // If not, or if node_n is hanging,
          // switch to the next neighbor.
          if (node_n == 8 || neighbor->is_hanging (node_n))
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
            case 4:
              if (node_side == 6)
                du_y.push_back (nedelec_gradient (neighbor, u, 4));
              break;
            case 5:
              if (node_side == 7)
                du_y.push_back (nedelec_gradient (neighbor, u, 5));
              break;
            case 6:
              if (node_side == 4)
                du_y.push_back (nedelec_gradient (neighbor, u, 4));
              break;
            case 7:
              if (node_side == 5)
                du_y.push_back (nedelec_gradient (neighbor, u, 5));
              break;  
            }

          if (weights_y.size () < du_y.size ())
            {
              weights_y[0] += 2 / hy;
              weights_y.push_back (-1 / hy);
            }
        }
    }

  //CCI
  // If on any boundary/interface...NON SO COME CHIAMARLA, AL MOMENTO
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

          //CCI
          // Check if neighbor contains the opposite vertex
          // of the "BOH" side containing "node".
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

          for (node_n = 0; node_n < 8; ++node_n)
            if (neighbor->gt (node_n) ==
                quadrant->gt (node_side))
              break;

          // If not, or if node_n is hanging,
          // switch to the next neighbor.
          if (node_n == 8 || neighbor->is_hanging (node_n))
            continue;

          hz = neighbor->p (2, 4) - neighbor->p (2, 0);

          switch (node_n)
            {
            case 0:
              if (node_side == 4)
                du_z.push_back (nedelec_gradient (neighbor, u, 8));
              break;
            case 1:
              if (node_side == 5)
                du_z.push_back (nedelec_gradient (neighbor, u, 9));
              break;
            case 2:
              if (node_side == 6)
                du_z.push_back (nedelec_gradient (neighbor, u, 10));
              break;
            case 3:
              if (node_side == 7)
                du_z.push_back (nedelec_gradient (neighbor, u, 11));
              break;
            case 4:
              if (node_side == 0)
                du_z.push_back (nedelec_gradient (neighbor, u, 8));
              break;
            case 5:
              if (node_side == 1)
                du_z.push_back (nedelec_gradient (neighbor, u, 9));
              break;
            case 6:
              if (node_side == 2)
                du_z.push_back (nedelec_gradient (neighbor, u, 10));
              break;
            case 7:
              if (node_side == 3)
                du_z.push_back (nedelec_gradient (neighbor, u, 11));
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

gradient3
bim2c_quadtree_pde_recovered_gradient (tmesh_3d& mesh,
                                       const q1_vec& u,
                                       active_fun3 is_active)
{
  q1_vec du_x_star (mesh.num_global_nodes (), 0);
  q1_vec du_y_star (mesh.num_global_nodes (), 0);
  q1_vec du_z_star (mesh.num_global_nodes (), 0);

  std::vector<bool> assigned_x (mesh.num_global_nodes (), 0);
  std::vector<bool> assigned_y (mesh.num_global_nodes (), 0);
  std::vector<bool> assigned_z (mesh.num_global_nodes (), 0);

  std::tuple<double, double, double,
         bool, bool, bool> du_star_loc; 

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
            bim2c_recovered_gradient_loc (quadrant, node,
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
/*
q2_vec3
bim2c_quadtree_pde_recovered_solution (tmesh_3d& mesh,
                                       const q1_vec& u,
                                       const gradient3& du)
{
  q2_vec3 u_star (mesh.num_local_quadrants (),
                 std::array<double, 27> ({0,0,0,0,0,0,0,0,0,
                 						  0,0,0,0,0,0,0,0,0,
                 						  0,0,0,0,0,0,0,0,0}));

  double hx = 0, hy = 0, hz = 0;

  std::array<double, 8> u_star_loc,
    					du_x_star_loc,
    					du_y_star_loc,
    					du_z_star_loc;

  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p (0, 1) - quadrant->p (0, 0);
      hy = quadrant->p (1, 2) - quadrant->p (1, 0);
      hz = quadrant->p (2, 4) - quadrant->p (2, 0);

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
        	  for (int pp = 0; pp < np; ++pp)
        	  {
        	  	u_star_loc[n] += u[quadrant->gparent(pp,n)];
        	  	du_x_star_loc[n] += std::get<0>(du)[quadrant->gparent (pp, n)];
        	  	du_y_star_loc[n] += std::get<1>(du)[quadrant->gparent (pp, n)];
        	  	du_z_star_loc[n] += std::get<2>(du)[quadrant->gparent (pp, n)];
        	  }
        	  u_star_loc[n] /= np;	
        	  du_x_star_loc[n] /= np;
        	  du_y_star_loc[n] /= np;
        	  du_z_star_loc[n] /= np;

        	  // Determine whether n is hanging on an edge
              // directed along the x or y or z direction.
              int i = 0; int p = 0;
              bool found = false;

              for (; i < 8 && !found; ++i)
                {
                  for (int pp = 0; pp < np && !found; ++pp)
                  	{
                  	  if (quadrant->parent (pp, n) == quadrant->t (i))	
                  	  	{
                  	  	  p = pp;
                  	  	  found = true;	
                  	  	}  	
                  	}
                }

              // Compute recovered solution at the
              // double-sized neighbor element.
              if (n == 0)
                {
                  if (i == 1)
                    u_star_loc[n] +=
                      (2 * hx) * (std::get<0>(du)[quadrant->gparent (, n)] -
                                  std::get<0>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 2)
                    u_star_loc[n] +=
                      (2 * hy) * (std::get<1>(du)[quadrant->gparent (, n)] -
                                  std::get<1>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 4)
                    u_star_loc[n] +=
                      (2 * hz) * (std::get<2>(du)[quadrant->gparent (, n)] -
                                  std::get<2>(du)[quadrant->gparent (, n)]) / 8;
                }  
              else if (n == 1)
                {
                  if (i == 0)
                    u_star_loc[n] +=
                      (2 * hx) * (std::get<0>(du)[quadrant->gparent (, n)] -
                                  std::get<0>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 3)
                    u_star_loc[n] +=
                      (2 * hy) * (std::get<1>(du)[quadrant->gparent (, n)] -
                                  std::get<1>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 5)
                    u_star_loc[n] +=
                      (2 * hz) * (std::get<2>(du)[quadrant->gparent (, n)] -
                                  std::get<2>(du)[quadrant->gparent (, n)]) / 8;
                } 
              else if (n == 2)
                {
                  if (i == 3)
                    u_star_loc[n] +=
                      (2 * hx) * (std::get<0>(du)[quadrant->gparent (, n)] -
                                  std::get<0>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 0)
                    u_star_loc[n] +=
                      (2 * hy) * (std::get<1>(du)[quadrant->gparent (, n)] -
                                  std::get<1>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 6)
                    u_star_loc[n] +=
                      (2 * hz) * (std::get<2>(du)[quadrant->gparent (, n)] -
                                  std::get<2>(du)[quadrant->gparent (, n)]) / 8;
                }
              else if (n == 3)
                {
                  if (i == 2)
                    u_star_loc[n] +=
                      (2 * hx) * (std::get<0>(du)[quadrant->gparent (, n)] -
                                  std::get<0>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 1)
                    u_star_loc[n] +=
                      (2 * hy) * (std::get<1>(du)[quadrant->gparent (, n)] -
                                  std::get<1>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 7)
                    u_star_loc[n] +=
                      (2 * hz) * (std::get<2>(du)[quadrant->gparent (, n)] -
                                  std::get<2>(du)[quadrant->gparent (, n)]) / 8;
                }
              else if (n == 4)
                {
                  if (i == 5)
                    u_star_loc[n] +=
                      (2 * hx) * (std::get<0>(du)[quadrant->gparent (, n)] -
                                  std::get<0>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 6)
                    u_star_loc[n] +=
                      (2 * hy) * (std::get<1>(du)[quadrant->gparent (, n)] -
                                  std::get<1>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 0)
                    u_star_loc[n] +=
                      (2 * hz) * (std::get<2>(du)[quadrant->gparent (, n)] -
                                  std::get<2>(du)[quadrant->gparent (, n)]) / 8;
                }          
              else if (n == 5)
                {
                  if (i == 4)
                    u_star_loc[n] +=
                      (2 * hx) * (std::get<0>(du)[quadrant->gparent (, n)] -
                                  std::get<0>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 7)
                    u_star_loc[n] +=
                      (2 * hy) * (std::get<1>(du)[quadrant->gparent (, n)] -
                                  std::get<1>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 1)
                    u_star_loc[n] +=
                      (2 * hz) * (std::get<2>(du)[quadrant->gparent (, n)] -
                                  std::get<2>(du)[quadrant->gparent (, n)]) / 8;
                }  
              else if (n == 6)
                {
                  if (i == 7)
                    u_star_loc[n] +=
                      (2 * hx) * (std::get<0>(du)[quadrant->gparent (, n)] -
                                  std::get<0>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 4)
                    u_star_loc[n] +=
                      (2 * hy) * (std::get<1>(du)[quadrant->gparent (, n)] -
                                  std::get<1>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 2)
                    u_star_loc[n] +=
                      (2 * hz) * (std::get<2>(du)[quadrant->gparent (, n)] -
                                  std::get<2>(du)[quadrant->gparent (, n)]) / 8;
                }    
              else if (n == 7)
                {
                  if (i == 6)
                    u_star_loc[n] +=
                      (2 * hx) * (std::get<0>(du)[quadrant->gparent (, n)] -
                                  std::get<0>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 5)
                    u_star_loc[n] +=
                      (2 * hy) * (std::get<1>(du)[quadrant->gparent (, n)] -
                                  std::get<1>(du)[quadrant->gparent (, n)]) / 8;
                  else if (i == 3)
                    u_star_loc[n] +=
                      (2 * hz) * (std::get<2>(du)[quadrant->gparent (, n)] -
                                  std::get<2>(du)[quadrant->gparent (, n)]) / 8;
                }   
          	}  

          u_star[quadrant->get_forest_quad_idx ()][n] = u_star_loc[n];  
        }
    }

  return u_star;  
}
*/
/* CCI: END ADDED */

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

/* CCI: BEGIN ADDED */

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

/*
  double result = (d11 * (y[1] - Y) * (x[1] - X) +
             d12 * (Y - y[0]) * (x[1] - X) +
             d21 * (y[1] - Y) * (X - x[0]) +
             d22 * (Y - y[0]) * (X - x[0])) / hxhyhz;
*/
  
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

  double num = u[0] * Xx1  * Yy1 * -Zz1 + 
           u[1] * -Xx0 * Yy1 * -Zz1 +
               u[2] * -Xx1 * Yy0 * -Zz1 + 
               u[3] * Xx0  * Yy0 * -Zz1 +
               u[4] * Xx1  * Yy1 * Zz0 + 
               u[5] * -Xx0 * Yy1 * Zz0 +
               u[6] * -Xx1 * Yy0 * Zz0 +
               u[7] * Xx0  * Yy0 * Zz0;          

  return (num / hxhyhz);
}

// Compute ||grad^* u - grad u||_L^2(q).
double
estimator_grad (tmesh_3d::quadrant_iterator q,
                const gradient3 & du_star,
                const q1_vec & u)
{
  double
    x[2] = {q->p(0,0), q->p(0,1)},
    y[2] = {q->p(1,0), q->p(1,3)},
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
    [x, y, z, dudxstar_loc, dudystar_loc, dudzstar_loc, u_loc]
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

/* CCI: END ADDED */