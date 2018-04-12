#include "quad_operators_3d.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <set>
#include <limits>
#include <iomanip>

void
bim3a_structure (tmesh &tmsh,
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

void
bim3a_dirichlet_bc (tmesh_3d& mesh, const dirichlet_bcs& bcs,
                    sparse_matrix& A, std::vector<double>& rhs)
{
  std::vector<double> row_sum (A.size ());
  
  // Set zero diagonal entries to sum (abs (row)).
  for (unsigned int row = 0; row < A.size (); ++row)
    {
      if (std::abs (A[row][row])
          < std::numeric_limits<double>::epsilon ())
        {
          row_sum[row] = std::accumulate
            (A[row].begin (), A[row].end (), 0.0,
              [] (double value,
                  const std::map<int, double>::value_type & p)
                {
                  return (value + std::abs (p.second));
                }
            );
        }
    }
  
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
                              rhs[col] -= A[col][row] * rhs[row];
                              A[col][row] = 0.0;
                            }
                        }
                    
                    if (std::abs (A[row][row])
                        < std::numeric_limits<double>::epsilon ())
                      A[row][row] = row_sum[row];
                    
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

