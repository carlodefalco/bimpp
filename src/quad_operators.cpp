#include "quad_operators.h"

#include <set>

void 
bim2a_advection_diffusion(tmesh & mesh,
                          const std::vector<double> & alpha,
                          const std::vector<double> & psi,
                          sparse_matrix & A)
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
            rows.push_back (quadrant->t(i));
          else
            {
              rows.push_back (quadrant->parent(0, i));
              rows.push_back (quadrant->parent(1, i));
            }
          
          for(int j = 0; j < 4; ++j)
            {
              cols.clear();
              
              if (!quadrant->is_hanging(j))
                cols.push_back (quadrant->t(j));
              else
                {
                  cols.push_back (quadrant->parent(0, j));
                  cols.push_back (quadrant->parent(1, j));
                }
              
              for (int r = 0; r < rows.size(); ++r)
                for (int c = 0; c < cols.size(); ++c)
                  {
                    A[rows[r]][cols[c]] += Aloc[i][j] / (rows.size() * cols.size());
                  }
            }
        }
    }
}

void bim2a_reaction (tmesh & mesh,
                     const std::vector<double> & delta,
                     const std::vector<double> & zeta,
                     sparse_matrix & A)
{
  double hx = 0, hy = 0;
  
  unsigned int iel = 0;
  std::vector<unsigned int> rows;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p(0, 1) - quadrant->p(0, 0);
      hy = quadrant->p(1, 2) - quadrant->p(1, 0);
      
      iel = quadrant->get_forest_quad_idx();
      
      for(int i = 0; i < 4; ++i)
        {
          rows.clear();
          
          if (!quadrant->is_hanging(i))
            rows.push_back (quadrant->t(i));
          else
            {
              rows.push_back (quadrant->parent(0, i));
              rows.push_back (quadrant->parent(1, i));
            }
          
          for (int r = 0; r < rows.size(); ++r)
            A[rows[r]][rows[r]] += (delta[iel] * zeta[quadrant->t(i)] * hx * hy / 4) / rows.size();
        }
    }
}

void bim2a_rhs (tmesh & mesh,
                const std::vector<double> & f,
                const std::vector<double> & g,
                std::vector<double> & rhs)
{
   double hx = 0, hy = 0;
   
   unsigned int iel = 0;
   std::vector<unsigned int> rows;
   
   for (auto quadrant = mesh.begin_quadrant_sweep ();
        quadrant != mesh.end_quadrant_sweep ();
        ++quadrant)
     {
        hx = quadrant->p(0, 1) - quadrant->p(0, 0);
        hy = quadrant->p(1, 2) - quadrant->p(1, 0);
        
        iel = quadrant->get_forest_quad_idx();
        
        for(int i = 0; i < 4; ++i)
          {
            rows.clear();
            
            if (!quadrant->is_hanging(i))
              rows.push_back (quadrant->t(i));
            else
              {
                rows.push_back (quadrant->parent(0, i));
                rows.push_back (quadrant->parent(1, i));
              }
            
            for (int r = 0; r < rows.size(); ++r)
              rhs[rows[r]] += (f[iel] * g[quadrant->t(i)] * hx * hy / 4) / rows.size();
          }
     }
}

void bim2a_dirichlet_bc (tmesh & mesh, const dirichlet_bcs & bcs,
                         sparse_matrix & A, std::vector<double> & rhs)
{
  int boundary_idx, tree_idx;
  unsigned int row, col;
  
  std::set<unsigned int> marked;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      tree_idx = quadrant->get_tree_idx();
      
      for (int i = 0; i < 4; ++i)
        {
          boundary_idx = quadrant->e(i);
          row = quadrant->t(i);
          
          // If current node is on boundary and has not been handled before.
          if (boundary_idx != tmesh::quadrant_t::NOT_ON_BOUNDARY
              && marked.count(row) == 0)
            {
              marked.insert(row); // Mark current node so to avoid duplicate operations.
              
              // Loop over all the boundary conditions.
              for (size_t bc = 0; bc < bcs.size(); ++bc)
                {
                  // If this boundary condition matches with the current node.
                  if (std::get<0>(bcs[bc]) == tree_idx
                      && std::get<1>(bcs[bc]) == boundary_idx)
                    {
                      // Impose boundary condition at rhs by evaluating it at the current node.
                      rhs[row] = (std::get<2>(bcs[bc]))(quadrant->p(0, i), quadrant->p(1, i));
                      
                      // Move non-diagonal entries from column "row" to rhs.
                      if (A[row].size ())
                        {
                          for (sparse_matrix::col_iterator j = A[row].begin ();
                               j != A[row].end ();
                               ++j)
                            {
                              col = A.col_idx(j);
                              
                              if (row != col)
                                {
                                  A[row][col] = 0.0;
                                  rhs[col] -= A[col][row] * rhs[row];
                                  A[col][row] = 0.0;
                                }
                            }
                        }
                      
                      // Multiply rhs by the diagonal entry.
                      rhs[row] *= A[row][row];
                    }
                }
            }
        }
    }
}
