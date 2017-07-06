#include <quad_operators.h>

void 
bim2a_advection_diffusion(tmesh & mesh,
                          const std::vector<double> & alpha,
                          const std::vector<double> & psi,
                          sparse_matrix & A)
{
  double psi12 = 0;
  double psi23 = 0;
  double psi43 = 0;
  double psi14 = 0;
  
  double bp12 = 0, bm12 = 0;
  double bp23 = 0, bm23 = 0;
  double bp43 = 0, bm43 = 0;
  double bp14 = 0, bm14 = 0;
  
  double hx = 0, hy = 0;
  
  std::array<std::array<double, 4>, 4> Aloc;
  
  unsigned int row = 0, col = 0;
  
  for (auto quadrant = mesh.begin_quadrant_sweep ();
       quadrant != mesh.end_quadrant_sweep ();
       ++quadrant)
    {
      psi12 = psi[quadrant->t(1)] - psi[quadrant->t(0)];
      psi23 = psi[quadrant->t(2)] - psi[quadrant->t(1)];
      psi43 = psi[quadrant->t(2)] - psi[quadrant->t(3)];
      psi14 = psi[quadrant->t(3)] - psi[quadrant->t(0)];
      
      bimu_bernoulli(psi12, bp12, bm12);
      bimu_bernoulli(psi23, bp23, bm23);
      bimu_bernoulli(psi43, bp43, bm43);
      bimu_bernoulli(psi14, bp14, bm14);
      
      hx = quadrant->p(0, 1) - quadrant->p(0, 0);
      hy = quadrant->p(1, 2) - quadrant->p(1, 0);
      
      bp12 *= alpha[quadrant->idx()] * hy / (4 * hx);
      bm12 *= alpha[quadrant->idx()] * hy / (4 * hx);
      bp23 *= alpha[quadrant->idx()] * hx / (4 * hy);
      bm23 *= alpha[quadrant->idx()] * hx / (4 * hy);
      bp43 *= alpha[quadrant->idx()] * hy / (4 * hx);
      bm43 *= alpha[quadrant->idx()] * hy / (4 * hx);
      bp14 *= alpha[quadrant->idx()] * hx / (4 * hy);
      bm14 *= alpha[quadrant->idx()] * hx / (4 * hy);
      
      Aloc[0][0] = bm12 + bm14;
      Aloc[0][1] = -bp12;
      Aloc[0][3] = -bp14;
      
      Aloc[1][0] = -bm12;
      Aloc[1][1] = bp12 + bm23;
      Aloc[1][2] = -bp23;
      
      Aloc[2][1] = -bm23;
      Aloc[2][2] = bp23 + bp43;
      Aloc[2][3] = -bm43;
        
      Aloc[3][0] = -bm14;
      Aloc[3][2] = -bp43;
      Aloc[3][3] = bm43 + bp14;
      
      for(int ii = 0; ii < 4; ++ii)
        {
          for(int jj = 0; jj < 4; ++jj)
            {
              row = quadrant->t(ii);
              col = quadrant->t(jj);
              
              A[row][col] += Aloc[ii][jj];
            }
        }
    }
}
