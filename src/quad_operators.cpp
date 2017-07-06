#include <quad_operators.h>

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
  
  std::array<std::array<double, 4>, 4> Aloc{};
  
  unsigned int iel = 0;
  unsigned int row = 0, col = 0;
  
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
      
      bp01 *= alpha[iel] * hy / (4 * hx);
      bm01 *= alpha[iel] * hy / (4 * hx);
      bp13 *= alpha[iel] * hx / (4 * hy);
      bm13 *= alpha[iel] * hx / (4 * hy);
      bp32 *= alpha[iel] * hy / (4 * hx);
      bm32 *= alpha[iel] * hy / (4 * hx);
      bp20 *= alpha[iel] * hx / (4 * hy);
      bm20 *= alpha[iel] * hx / (4 * hy);
      
      Aloc[0][0] = bm01 + bp20;
      Aloc[0][1] = -bp01;
      Aloc[0][2] = -bm20;
      
      Aloc[1][0] = -bm01;
      Aloc[1][1] = bp01 + bm13;
      Aloc[1][3] = -bp13;
      
      Aloc[2][0] = -bp20;
      Aloc[2][2] = bp32 + bm20;
      Aloc[2][3] = -bm32;
      
      Aloc[3][1] = -bm13;
      Aloc[3][2] = -bp32;
      Aloc[3][3] = bm32 + bp13;
      
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

void bim2a_rhs (tmesh & mesh,
                const std::vector<double> & f,
                const std::vector<double> & g,
                std::vector<double> & rhs)
{
   double hx = 0, hy = 0;
   
   unsigned int iel = 0;
   unsigned int row = 0;
   
   for (auto quadrant = mesh.begin_quadrant_sweep ();
        quadrant != mesh.end_quadrant_sweep ();
        ++quadrant)
     {
        hx = quadrant->p(0, 1) - quadrant->p(0, 0);
        hy = quadrant->p(1, 2) - quadrant->p(1, 0);
        
        iel = quadrant->get_forest_quad_idx();
        
        for(int ii = 0; ii < 4; ++ii)
        {
           row = quadrant->t(ii);
           
           rhs[row] += f[iel] * g[quadrant->t(ii)] * hx * hy / 4;
        }
     }
}
