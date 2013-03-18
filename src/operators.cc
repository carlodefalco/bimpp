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

void bim3a_structure (const mesh& msh, sparse_matrix& SG)
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

void bim3a_rhs       (const mesh& msh, const std::vector<double>& ecoeff, const std::vector<double>& ncoeff, std::vector<double>& b)
{
  b.resize (msh.nnodes);
  for (int iel = 0; iel < msh.nelements; ++iel)
    for (int inode = 0; inode < 4; ++inode)
      {
        int ig = msh.t (inode, iel);
        b[ig] += ncoeff[ig] * (ecoeff[iel] / 4.0) * msh.volume (iel);
      }
};

void bim3a_reaction  (const mesh& msh, const std::vector<double>& ecoeff, const std::vector<double>& ncoeff, sparse_matrix& A)
{
  if (A.size () < msh.nnodes)
    A.resize (msh.nnodes);
  
  int iel, inode, jnode;

  for (iel = 0; iel < msh.nelements; ++iel)
    for (inode = 0; inode < 4; ++inode)
      {
        int ig = msh.t (inode, iel);
        
        if (ig >= msh.nnodes)
          std::cout << " out of bounds" << std::endl;
        
        A[ig][ig] += ncoeff[ig] * ecoeff[iel] * msh.volume (iel) / 4.0;
      }
};

void bim3a_laplacian (const mesh& msh, const std::vector<double>& acoeff, sparse_matrix& SG)
{

  if (SG.size () < msh.nnodes)
    SG.resize (msh.nnodes);

  double epsilonareak;  
  int iel, inode, jnode, idir;

  for (iel = 0; iel < msh.nelements; ++iel)
    {      
      epsilonareak = acoeff[iel] * msh.volume (iel);
  
      // Compute local laplacian matrix and assemble into global matrix
      for (inode = 0; inode < 4; ++inode)
        for (jnode = 0; jnode < 4; ++jnode)
          {
            int ig = msh.t (inode, iel);
            int jg = msh.t (jnode, iel);

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

void bim3a_advection_diffusion (const mesh& msh, const std::vector<double>& acoeff, const std::vector<double>& v, sparse_matrix& SG)
{

  if (SG.size () < msh.nnodes)
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

            if (ginode[inode][jnode] >= msh.nnodes || gjnode[inode][jnode] >= msh.nnodes)
              std::cout << " out of bounds" << std::endl;

            for (idir = 0; idir < 3; ++idir)
              Lloc[inode][jnode] += msh.shg (idir, inode, iel) * msh.shg (idir, jnode, iel) * epsilonareak;
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


void bimu_bernoulli (double x, double &bp, double &bn)
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
      while (abs (df) > 1.0e-16)
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
