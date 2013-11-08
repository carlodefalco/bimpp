// Copyright (C) 2004-2011  Carlo de Falco
//
// This file is part of:
//     secs3d - A 3-D Drift--Diffusion Semiconductor Device Simulator 
//
//  secs3d is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  secs3d is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with secs3d; If not, see <http://www.gnu.org/licenses/>.
//
//  author: Carlo de Falco     <cdf _AT_ users.sourceforge.net>

#include <operators.h>
#include <cstring>

//namespace bim
//{

void 
bim3a_structure (const mesh& msh, sparse_matrix& SG)
{
  SG.resize (msh.nnodes);
  for (int iel = 0; iel < msh.nelements; ++iel)
    for (int inode = 0; inode < 4; ++inode)
      for (int jnode = 0; jnode < 4; ++jnode)
        {
          int ig = msh.t (inode, iel);
          int jg = msh.t (jnode, iel);
          SG[ig][jg] = 0.0;
        }
  SG.set_properties ();
};

void 
bim3a_rhs  
(const mesh& msh, const std::vector<double>& ecoeff, 
 const std::vector<double>& ncoeff, std::vector<double>& b)
{
  b.resize (msh.nnodes);
  for (int iel = 0; iel < msh.nelements; ++iel)
    for (int inode = 0; inode < 4; ++inode)
      {
        int ig = msh.t (inode, iel);
        b[ig] += ncoeff[ig] * (ecoeff[iel] / 4.0) * msh.volume (iel);
      }
};

void 
bim3a_reaction  
(const mesh& msh, const std::vector<double>& ecoeff, 
 const std::vector<double>& ncoeff, sparse_matrix& A)
{
  if (A.size () < size_t (msh.nnodes))
    A.resize (msh.nnodes);
  
  int iel, inode;

  for (iel = 0; iel < msh.nelements; ++iel)
    for (inode = 0; inode < 4; ++inode)
      {
        int ig = msh.t (inode, iel);
        
        if (ig >= msh.nnodes)
          std::cout << " out of bounds" << std::endl;
        
        A[ig][ig] += ncoeff[ig] * ecoeff[iel] * msh.volume (iel) / 4.0;
      }
};
  
void
bim3a_local_laplacian (const mesh& msh, 
                       const int iel, 
                       const double acoeff, 
                       double SG[16])
{
  int inode, jnode, idir;
  double ishg, jshg;
#define sgloc(i, j) (SG[i + 4*j])
    
  for (inode = 0; inode < 4; ++inode)
    for (jnode = 0; jnode < 4; ++jnode)
      {
        for (idir = 0; idir < 3; ++idir)
          {
            ishg = msh.shg (idir, inode, iel);
            jshg = msh.shg (idir, jnode, iel);
            sgloc(inode, jnode) += ishg * jshg *  acoeff * msh.volume (iel);
          }
      }
#undef sgloc
};

void 
bim3a_laplacian (const mesh& msh, 
                 const std::vector<double>& acoeff, 
                 sparse_matrix& SG)
{
    
  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double epsilonareak;  
  int iel, inode, jnode, 
    idir, ig, jg;

  for (iel = 0; iel < msh.nelements; ++iel)
    {      
      epsilonareak = acoeff[iel] * msh.volume (iel);
  
      // Compute local laplacian matrix and assemble into global matrix
      for (inode = 0; inode < 4; ++inode)
        for (jnode = 0; jnode < 4; ++jnode)
          {
            ig = msh.t (inode, iel);
            jg = msh.t (jnode, iel);

            if (ig >= msh.nnodes || jg >= msh.nnodes)
              std::cout << " out of bounds" << std::endl;

            for (idir = 0; idir < 3; ++idir)
              {
                double ishg = msh.shg (idir, inode, iel);
                double jshg = msh.shg (idir, jnode, iel);
                SG[ig][jg] += ishg * jshg * epsilonareak;
              }
          }
    }
};

void
bim3a_advection_diffusion (const mesh& msh, 
                           const std::vector<double>& acoeff, 
                           const std::vector<double>& v, 
                           sparse_matrix& SG)
{

  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);
    
  double Lloc[4][4], Sloc[4][4];
  double epsilonareak;
  double 
    bm12, bm13, bm14, bm23, bm24, bm34,
    bp12, bp13, bp14, bp23, bp24, bp34;

  int ginode[4][4], gjnode[4][4];
  double vloc[4];

  int iel, inode, jnode, idir, ig, jg;

  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (inode = 0; inode < 4; ++inode)
        for (jnode = 0; jnode < 4; ++jnode)
          Lloc[inode][jnode] = 0.0;
      
      epsilonareak = acoeff[iel] * msh.volume (iel);
  
      // Compute local laplacian matrix
      for (inode = 0; inode < 4; ++inode)
        for (jnode = 0; jnode < 4; ++jnode)
          {
            ginode[inode][jnode] = msh.t (inode, iel);
            gjnode[inode][jnode] = msh.t (jnode, iel);

            if (ginode[inode][jnode] >= msh.nnodes 
                || gjnode[inode][jnode] >= msh.nnodes)
              std::cout << " out of bounds" << std::endl;

            for (idir = 0; idir < 3; ++idir)
              Lloc[inode][jnode] += msh.shg (idir, inode, iel) * 
                msh.shg (idir, jnode, iel) * epsilonareak;
          }


      for (inode = 0; inode < 4; ++inode)
        vloc[inode] = v[msh.t (inode, iel)];
        
      bimu_bernoulli (vloc[1]-vloc[0], bp12, bm12);
      bimu_bernoulli (vloc[2]-vloc[0], bp13, bm13);
      bimu_bernoulli (vloc[3]-vloc[0], bp14, bm14);
      bimu_bernoulli (vloc[2]-vloc[1], bp23, bm23);
      bimu_bernoulli (vloc[3]-vloc[1], bp24, bm24);
      bimu_bernoulli (vloc[3]-vloc[2], bp34, bm34);

      bp12 *= Lloc[0][1];
      bm12 *= Lloc[0][1];
      bp13 *= Lloc[0][2];
      bm13 *= Lloc[0][2];
      bp14 *= Lloc[0][3];
      bm14 *= Lloc[0][3];
      bp23 *= Lloc[1][2];
      bm23 *= Lloc[1][2];
      bp24 *= Lloc[1][3];
      bm24 *= Lloc[1][3];
      bp34 *= Lloc[2][3];
      bm34 *= Lloc[2][3];
        
      /*
        ## Sloc=[...
        ##        -bm12-bm13-bm14,bp12            ,bp13           ,bp14     
        ##        bm12           ,-bp12-bm23-bm24 ,bp23           ,bp24
        ##        bm13           ,bm23            ,-bp13-bp23-bm34,bp34
        ##        bm14           ,bm24            ,bm34           ,-bp14-bp24-bp34...
        ##       ];
      */
      
      Sloc[0][0] = -bm12-bm13-bm14;
      Sloc[0][1] = bp12;
      Sloc[0][2] = bp13;
      Sloc[0][3] = bp14;

      Sloc[1][0] = bm12;
      Sloc[1][1] = -bp12-bm23-bm24; 
      Sloc[1][2] = bp23;
      Sloc[1][3] = bp24;

      Sloc[2][0] = bm13;
      Sloc[2][1] = bm23;
      Sloc[2][2] = -bp13-bp23-bm34;
      Sloc[2][3] = bp34;
  
      Sloc[3][0] = bm14;
      Sloc[3][1] = bm24;
      Sloc[3][2] = bm34;
      Sloc[3][3] = -bp14-bp24-bp34;

      // assemble global matrix
      for (inode = 0; inode < 4; ++inode)
        for (jnode = 0; jnode < 4; ++jnode)
          {
            ig = ginode[inode][jnode];
            jg = gjnode[inode][jnode];
                     
            SG[ig][jg] += Sloc[inode][jnode];            
          }
    }
};


void
bim3a_osc_local_laplacian (const mesh& msh, const int iel, 
                           const double acoeff, double *Lloc)
{
  int inode, idir;
  double A[12] = {0}, Ann[4]= {0}, AidotAj[6]={0}, r[12] = {0};
  double vol = msh.volume (iel);
  double epsilonareak  = acoeff / vol / 48.0;  

  for (inode = 0; inode < 4; ++inode)
    {
      A[0 + 4 * inode] = 3.0 * vol * msh.shg (0, inode, iel);
      A[1 + 4 * inode] = 3.0 * vol * msh.shg (1, inode, iel);
      A[2 + 4 * inode] = 3.0 * vol * msh.shg (2, inode, iel);

      Ann[inode] = pow (A[0 + 4 * inode], 2) +
        pow (A[1 + 4 * inode], 2) + 
        pow (A[2 + 4 * inode], 2);

      memset (&(r[3*inode]), 0, 3*sizeof(double));
    }
  memset (&(AidotAj[0]), 0, 6 * sizeof(double));

  for (idir = 0; idir < 3; ++idir) 
    {
      r[0]  += (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (0 , iel))) * 
        (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (1 , iel))); //rik dot rjk 
      r[1]  += (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (0 , iel))) * 
        (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (1 , iel))); //ril dot rjl 
      r[2]  += (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (0 , iel))) * 
        (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (2 , iel))); //rij dot rkj 
      r[3]  += (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (0 , iel))) * 
        (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (2 , iel))); //ril dot rkl 
      r[4]  += (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (0 , iel))) * 
        (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (3 , iel))); //rij dot rlj 
      r[5]  += (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (0 , iel))) * 
        (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (3 , iel))); //rik dot rlk 
      r[6]  += (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (1 , iel))) * 
        (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (2 , iel))); //rji dot rki 
      r[7]  += (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (1 , iel))) * 
        (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (2 , iel))); //rjl dot rkl 
      r[8]  += (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (1 , iel))) * 
        (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (3 , iel))); //rji dot rli 
      r[9]  += (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (1 , iel))) * 
        (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (3 , iel))); //rjk dot rlk 
      r[10] += (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (2 , iel))) * 
        (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (3 , iel))); //rki dot rli 
      r[11] += (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (2 , iel))) * 
        (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (3 , iel))); //rkj dot rlj 

      AidotAj[0] += A[idir + 4 * 2] * A[idir + 4 * 3]; // Ak dot Al
      AidotAj[1] += A[idir + 4 * 1] * A[idir + 4 * 3]; // Aj dot Al
      AidotAj[2] += A[idir + 4 * 1] * A[idir + 4 * 2]; // Aj dot Ak
      AidotAj[3] += A[idir + 4 * 0] * A[idir + 4 * 3]; // Ai dot Al
      AidotAj[4] += A[idir + 4 * 0] * A[idir + 4 * 2]; // Ai dot Ak
      AidotAj[5] += A[idir + 4 * 0] * A[idir + 4 * 1]; // Ai dot Aj
    }

  Lloc[0 + 4 * 1] = Lloc[1 + 4 * 0] = - epsilonareak * (2.0 * r[0]  * r[1]  + AidotAj[0] * (r[0]  * r[0]  / Ann[3] + r[1]  * r[1]  / Ann[2]));
  Lloc[0 + 4 * 2] = Lloc[2 + 4 * 0] = - epsilonareak * (2.0 * r[2]  * r[3]  + AidotAj[1] * (r[2]  * r[2]  / Ann[3] + r[3]  * r[3]  / Ann[1]));
  Lloc[0 + 4 * 3] = Lloc[3 + 4 * 0] = - epsilonareak * (2.0 * r[4]  * r[5]  + AidotAj[2] * (r[4]  * r[4]  / Ann[2] + r[5]  * r[5]  / Ann[1]));
  Lloc[1 + 4 * 2] = Lloc[2 + 4 * 1] = - epsilonareak * (2.0 * r[6]  * r[7]  + AidotAj[3] * (r[6]  * r[6]  / Ann[3] + r[7]  * r[7]  / Ann[0]));
  Lloc[1 + 4 * 3] = Lloc[3 + 4 * 1] = - epsilonareak * (2.0 * r[8]  * r[9]  + AidotAj[4] * (r[8]  * r[8]  / Ann[2] + r[9]  * r[9]  / Ann[0]));
  Lloc[2 + 4 * 3] = Lloc[3 + 4 * 2] = - epsilonareak * (2.0 * r[10] * r[11] + AidotAj[5] * (r[10] * r[10] / Ann[1] + r[11] * r[11] / Ann[0]));
  Lloc[0 + 4 * 0] = - Lloc[0 + 4 * 1] - Lloc[0 + 4 * 2] - Lloc[0 + 4 * 3];
  Lloc[1 + 4 * 1] = - Lloc[1 + 4 * 0] - Lloc[1 + 4 * 2] - Lloc[1 + 4 * 3];
  Lloc[2 + 4 * 2] = - Lloc[2 + 4 * 0] - Lloc[2 + 4 * 1] - Lloc[2 + 4 * 3];
  Lloc[3 + 4 * 3] = - Lloc[3 + 4 * 0] - Lloc[3 + 4 * 1] - Lloc[3 + 4 * 2];

};


void
bim3a_osc_laplacian (const mesh& msh, 
                     const std::vector<double>& acoeff, 
                     sparse_matrix& SG)
{

  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16];
  int iel, inode, jnode;

  for (iel = 0; iel < msh.nelements; ++iel)
    {      
       
      // Compute local laplacian matrix and assemble into global matrix
      memset (Lloc, 0, 16*sizeof(double));
      bim3a_osc_local_laplacian (msh,  iel, acoeff[iel], Lloc);

      for (inode = 0; inode < 4; ++inode)
        for (jnode = 0; jnode < 4; ++jnode)
          {
            int ig = msh.t (inode, iel);
            int jg = msh.t (jnode, iel);

            if (ig >= msh.nnodes || jg >= msh.nnodes)
              std::cout << " out of bounds" << std::endl;

            SG[ig][jg] += Lloc[inode + 4 * jnode];
          }
    }
      
};

void
bim3a_osc_local_advection_diffusion (const mesh& msh, const int iel, 
                                     const double acoeff, double *Sloc)
{    
  // int ginode[16], gjnode[16];
  // int iel, inode, jnode, idir, ig, jg;

  // double vloc[4] = {0}, 
  //   bm[6]= {0}, 
  //     bp[6]= {0}
  //     Lloc = {0};

  // double 
  //   Sloc [16] = {0},       
  //   A[12] = {0}, 
  //     Ann[4] = {0}, 
  //       AidotAj[6] = {0}, 
  //         r[12];

  //   double epsilonareak;  
};

void
bim3a_osc_advection_diffusion (const mesh& msh, 
                               const std::vector<double>& acoeff, 
                               const std::vector<double>& v, 
                               sparse_matrix& SG)
{

  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);
    
  int ginode[4][4], gjnode[4][4];
  int iel, inode, jnode, idir, ig, jg;

  double vloc[4], bm[6], bp[6], Lloc[4][4], Sloc[4][4], 
    A[3][4], Ann[4], AidotAj[6], r[12];
  double epsilonareak;  

  for (iel = 0; iel < msh.nelements; ++iel)
    {      
      epsilonareak = acoeff[iel] / msh.volume (iel) / 4.8e1;
        
      for (inode = 0; inode < 4; ++inode)
        {
          Lloc[inode][0] = 0;
          Lloc[inode][1] = 0;
          Lloc[inode][2] = 0;
          Lloc[inode][3] = 0;
          A[0][inode] = 3 * msh.volume (iel) * msh.shg (0, inode, iel);
          A[1][inode] = 3 * msh.volume (iel) * msh.shg (1, inode, iel);
          A[2][inode] = 3 * msh.volume (iel) * msh.shg (2, inode, iel);
          Ann[inode] = A[0][inode] * A[0][inode] + 
            A[1][inode] * A[1][inode] + 
            A[2][inode] * A[2][inode];
          r[3 * inode] = r[3 * inode + 1] = r[3 * inode + 2] = 0.0; 
          vloc[inode] = v[msh.t (inode, iel)];
        }
      AidotAj[0] = AidotAj[1] = AidotAj[2] = AidotAj[3] = AidotAj[4] = AidotAj[5] = 0.0; 
      //bm[0] = bm[1] = bm[2] = bm[3] = bm[4] = bm[5] = 0.0; 
      //bp[0] = bp[1] = bp[2] = bp[3] = bp[4] = bp[5] = 0.0; 

      for (idir = 0; idir < 3; ++idir) 
        {
          r[0]  += (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (0 , iel))) *
            (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (1 , iel))); //rik dot rjk
          r[1]  += (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (0 , iel))) *
            (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (1 , iel))); //ril dot rjl
          r[2]  += (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (0 , iel))) *
            (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (2 , iel))); //rij dot rkj
          r[3]  += (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (0 , iel))) *
            (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (2 , iel))); //ril dot rkl
          r[4]  += (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (0 , iel))) *
            (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (3 , iel))); //rij dot rlj
          r[5]  += (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (0 , iel))) *
            (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (3 , iel))); //rik dot rlk
          r[6]  += (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (1 , iel))) *
            (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (2 , iel))); //rji dot rki
          r[7]  += (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (1 , iel))) *
            (msh.p(idir, msh.t (3 , iel)) - msh.p(idir, msh.t (2 , iel))); //rjl dot rkl
          r[8]  += (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (1 , iel))) *
            (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (3 , iel))); //rji dot rli
          r[9]  += (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (1 , iel))) *
            (msh.p(idir, msh.t (2 , iel)) - msh.p(idir, msh.t (3 , iel))); //rjk dot rlk
          r[10] += (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (2 , iel))) *
            (msh.p(idir, msh.t (0 , iel)) - msh.p(idir, msh.t (3 , iel))); //rki dot rli
          r[11] += (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (2 , iel))) *
            (msh.p(idir, msh.t (1 , iel)) - msh.p(idir, msh.t (3 , iel))); //rkj dot rlj
          AidotAj[0] += A[idir][2] * A[idir][3]; // Ak dot Al
          AidotAj[1] += A[idir][1] * A[idir][3]; // Aj dot Al
          AidotAj[2] += A[idir][1] * A[idir][2]; // Aj dot Ak
          AidotAj[3] += A[idir][0] * A[idir][3]; // Ai dot Al
          AidotAj[4] += A[idir][0] * A[idir][2]; // Ai dot Ak
          AidotAj[5] += A[idir][0] * A[idir][1]; // Ai dot Aj
        }


      Lloc[0][1] =  - epsilonareak * (2.0 * r[0]  * r[1]  + AidotAj[0] * (r[0]  * r[0]  / Ann[3] + r[1]  * r[1]  / Ann[2]));
      Lloc[0][2] =  - epsilonareak * (2.0 * r[2]  * r[3]  + AidotAj[1] * (r[2]  * r[2]  / Ann[3] + r[3]  * r[3]  / Ann[1]));
      Lloc[0][3] =  - epsilonareak * (2.0 * r[4]  * r[5]  + AidotAj[2] * (r[4]  * r[4]  / Ann[2] + r[5]  * r[5]  / Ann[1]));
      Lloc[1][2] =  - epsilonareak * (2.0 * r[6]  * r[7]  + AidotAj[3] * (r[6]  * r[6]  / Ann[3] + r[7]  * r[7]  / Ann[0]));
      Lloc[1][3] =  - epsilonareak * (2.0 * r[8]  * r[9]  + AidotAj[4] * (r[8]  * r[8]  / Ann[2] + r[9]  * r[9]  / Ann[0]));
      Lloc[2][3] =  - epsilonareak * (2.0 * r[10] * r[11] + AidotAj[5] * (r[10] * r[10] / Ann[1] + r[11] * r[11] / Ann[0]));
      //Lloc[1][0] = Lloc[0][1];
      //Lloc[2][0] = Lloc[0][2];
      //Lloc[3][0] = Lloc[0][3];
      //Lloc[2][1] = Lloc[1][2];
      //Lloc[3][1] = Lloc[1][3];
      //Lloc[3][2] = Lloc[2][3];
      Lloc[0][0] = - Lloc[0][1] - Lloc[0][2] - Lloc[0][3];
      Lloc[1][1] = - Lloc[1][0] - Lloc[1][2] - Lloc[1][3];
      Lloc[2][2] = - Lloc[2][0] - Lloc[2][1] - Lloc[2][3];
      Lloc[3][3] = - Lloc[3][0] - Lloc[3][1] - Lloc[3][2];
  
      bimu_bernoulli (vloc[1]-vloc[0], bp[0], bm[0]);
      bimu_bernoulli (vloc[2]-vloc[0], bp[1], bm[1]);
      bimu_bernoulli (vloc[3]-vloc[0], bp[2], bm[2]);
      bimu_bernoulli (vloc[2]-vloc[1], bp[3], bm[3]);
      bimu_bernoulli (vloc[3]-vloc[1], bp[4], bm[4]);
      bimu_bernoulli (vloc[3]-vloc[2], bp[5], bm[5]);

      bp[0] *= Lloc[0][1];        bm[0] *= Lloc[0][1];
      bp[1] *= Lloc[0][2];        bm[1] *= Lloc[0][2];
      bp[2] *= Lloc[0][3];        bm[2] *= Lloc[0][3];
      bp[3] *= Lloc[1][2];        bm[3] *= Lloc[1][2];
      bp[4] *= Lloc[1][3];        bm[4] *= Lloc[1][3];
      bp[5] *= Lloc[2][3];        bm[5] *= Lloc[2][3];
        
      /*
        ## Sloc=[-bm0-bm1-bm2,bp0          ,bp1          ,bp2
        ##       bm0         ,-bp0-bm3-bm4 ,bp3          ,bp4
        ##       bm1         ,bm3          ,-bp1-bp3-bm5 ,bp5
        ##       bm2         ,bm4          ,bm5          ,-bp2-bp4-bp5 ];
      */
      
      Sloc[0][0] = -bm[0] -bm[1] -bm[2];
      Sloc[0][1] = bp[0];
      Sloc[0][2] = bp[1];
      Sloc[0][3] = bp[2];

      Sloc[1][0] = bm[0];
      Sloc[1][1] = -bp[0] -bm[3] -bm[4]; 
      Sloc[1][2] = bp[3];
      Sloc[1][3] = bp[4];

      Sloc[2][0] = bm[1];
      Sloc[2][1] = bm[3];
      Sloc[2][2] = -bp[1] -bp[3] -bm[5];
      Sloc[2][3] = bp[5];
  
      Sloc[3][0] = bm[2];
      Sloc[3][1] = bm[4];
      Sloc[3][2] = bm[5];
      Sloc[3][3] = -bp[2] -bp[4] -bp[5];

      // assemble global matrix
      for (inode = 0; inode < 4; ++inode)
        for (jnode = 0; jnode < 4; ++jnode)
          {
            ig = ginode[inode][jnode];
            jg = gjnode[inode][jnode];
                     
            SG[ig][jg] += Sloc[inode][jnode];            
          }
    }
      
  //std::cout << std::endl;
  //std::cout << " ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ " << std::endl;
  //std::cout << " bim3a_osc_advection_diffusion(...): " << std::endl;
  //std::cout << " TO DO ! " << std::endl;
  //std::cout << " ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ " << std::endl;
};

void 
bimu_bernoulli (double x, double &bp, double &bn)
{
  const double xlim= 1.0e-2;
  double ax  = fabs (x);

  bp  = 0.0;
  bn  = 0.0;
  
  //  X=0
  if (x == 0.0)
    {
      bp = 1.0;
      bn = 1.0;
      return;
    }
  
  // ASYMPTOTICS
  if (ax > 80.0)
    {
      if (x > 0.0)
        {
          bp = 0.0;
          bn = x;
        }
      else
        {
          bp = -x;
          bn = 0.0;
        }
      return;
    }
  
  // INTERMEDIATE VALUES
  if (ax <= 80 &&  ax > xlim)
    {
      bp = x / (exp (x) - 1.0);
      bn = x + bp;
      return;
    }

  // SMALL VALUES
  if (ax <= xlim &&  ax != 0.0)
    {
      double jj = 1.0;
      double fp = 1.0;
      double fn = 1.0;
      double df = 1.0;
      double segno = 1.0;
      while (fabs (df) > 1.0e-16)
        {
          jj += 1.0;
          segno = -segno;
          df = df * x / jj;
          fp = fp + df;
          fn = fn + segno * df;
        }
      bp = 1 / fp;
      bn = 1 / fn;
      return;
    }
 
};
//}
