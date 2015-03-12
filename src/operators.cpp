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

/*! \file operators.cpp
  \brief Functions to build the discrete version of differential operators.
*/
#include <operators.h>
#include <cstring>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
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
bim3a_rhs (mesh& msh,
           const std::vector<double>& ecoeff,
           const std::vector<double>& ncoeff,
           std::vector<double>& b)
{
  if (b.size () < size_t (msh.nnodes))
    b.resize (msh.nnodes);
  for (int iel = 0; iel < msh.nelements; ++iel)
    for (int inode = 0; inode < 4; ++inode)
      {
        int ig = msh.t (inode, iel);
        b[ig] += ncoeff[ig] * (ecoeff[iel] / 4.0) * msh.volume (iel);
      }
};

void
bim3a_reaction (mesh& msh,
                const std::vector<double>& ecoeff,
                const std::vector<double>& ncoeff,
                sparse_matrix& A)
{
  if (A.size () < size_t (msh.nnodes))
    A.resize (msh.nnodes);

  int iel, inode[4];
  double Lloc[16], nloc[4] = {0.0, 0.0, 0.0, 0.0};

  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          inode[ii] = msh.t (ii, iel);
          nloc[ii] = ncoeff[inode[ii]];
        }

      memset (Lloc, 0, 16 * sizeof (double));
      bim3a_local_reaction (& (msh.shp (0, 0)),
                            & (msh.wjacdet (0, iel)),
                            ecoeff[iel],
                            nloc, Lloc);

      for (int ii = 0; ii < 4; ++ii)
        A[inode[ii]][inode[ii]] += Lloc[5*ii];

    }
};

void
bim3a_laplacian (mesh& msh,
                 const std::vector<double>& acoeff,
                 sparse_matrix& SG)
{
  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16];
  int iel, inode[4];

  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        inode[ii] = msh.t (ii, iel);

      memset (Lloc, 0, 16 * sizeof (double));
      bim3a_local_laplacian (& (msh.shg (0, 0, iel)),
                             msh.volume (iel),
                             acoeff[iel],
                             Lloc);

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          SG[inode[ii]][inode[jj]] += Lloc[ii + 4 * jj];
    }
};

void
bim3a_laplacian_anisotropic (mesh& msh,
                             const std::vector<double>& acoeff,
                             sparse_matrix& SG)
{
  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16];
  int iel, inode[4];

  double adim = (double) acoeff.size () / msh.nelements;
  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        inode[ii] = msh.t (ii, iel);

      memset (Lloc, 0, 16 * sizeof (double));

      if (adim == 3)
        bim3a_local_laplacian_anisotropic_diag (& (msh.shg (0, 0, iel)),
                                                msh.volume (iel),
                                                &acoeff[iel*3],
                                                Lloc);
      else if (adim == 9)
        bim3a_local_laplacian_anisotropic (& (msh.shg (0, 0, iel)),
                                           msh.volume (iel),
                                          &acoeff[iel*9],
                                          Lloc);
      else
        {
          std::cerr << "bim3a_laplacian_anisotropic: "
                    << "coefficient acoeff has wrong dimension" << std::endl;
          exit (-1);
        }

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          SG[inode[ii]][inode[jj]] += Lloc[ii + 4 * jj];
    }
};

void
bim3a_advection_diffusion (mesh& msh,
                           const std::vector<double>& acoeff,
                           const std::vector<double>& v,
                           sparse_matrix& SG)
{
  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16], vloc[4];
  int iel, inode[4];

  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          inode[ii] = msh.t (ii, iel);
          vloc[ii] = v[inode[ii]];
        }

      memset (Lloc, 0, 16 * sizeof (double));
      bim3a_local_laplacian (& (msh.shg (0, 0, iel)),
                             msh.volume (iel),
                             acoeff[iel],
                             Lloc);

      bim3a_local_advection (vloc, 1.0, 1.0, Lloc);

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          SG[inode[ii]][inode[jj]] += Lloc[ii + 4 * jj];
    }
};

void
bim3a_advection_diffusion_anisotropic (mesh& msh,
                                       const std::vector<double>& acoeff,
                                       const std::vector<double>& v,
                                       sparse_matrix& SG)
{
  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16], vloc[4];
  int iel, inode[4];

  double adim = (double) acoeff.size () / msh.nelements;

  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          inode[ii] = msh.t (ii, iel);
          vloc[ii] = v[inode[ii]];
        }

      memset (Lloc, 0, 16 * sizeof (double));
      if (adim == 3)
        bim3a_local_laplacian_anisotropic_diag (& (msh.shg (0, 0, iel)),
                                                msh.volume (iel),
                                                &acoeff[iel*3],
                                                Lloc);

      else if (adim == 9)
        bim3a_local_laplacian_anisotropic (& (msh.shg (0, 0, iel)),
                                           msh.volume (iel),
                                           &acoeff[iel*9],
                                           Lloc);
      else
        {
          std::cerr << "bim3a_advection_diffusion_anisotropic: "
                    << "coefficient acoeff has wrong dimension" << std::endl;
          exit (-1);
        }

      bim3a_local_advection (vloc, 1.0, 1.0, Lloc);

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          SG[inode[ii]][inode[jj]] += Lloc[ii + 4 * jj];
    }
};


void
bim3a_advection_upwind (mesh &msh,
                        const std::vector<double>& v,
                        sparse_matrix& UP)
{
  if (UP.size () < size_t (msh.nnodes))
    UP.resize (msh.nnodes);

  double Lloc[16];
  int iel, inode[4];

  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          inode[ii] = msh.t (ii, iel);
        }

      memset (Lloc, 0, 16 * sizeof (double));
      bim3a_local_laplacian (& (msh.shg (0, 0, iel)),
                             msh.volume (iel),
                             1.0,
                             Lloc);

      double x[4], y[4], z[4];
      double vnodes[6]; //v12,v13,v14,v23,v24,v34;
      for (int i = 0; i < 4; ++i)
        {
          x[i] = msh.p (0, inode[i]);
          y[i] = msh.p (1, inode[i]);
          z[i] = msh.p (2, inode[i]);
        }

      if (v.size () == 1)
        for (int i = 0; i < 6; ++i)
          vnodes[i] = 0;

      else if (v.size () == 3u * msh.nelements)
        {
          int node = 0;
          for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 4; ++j){
              vnodes[node] = v[iel * 3] * (x[j] - x[i]) +
                             v[iel * 3 + 1] * (y[j] - y[i]) +
                             v[iel * 3 + 2] * (z[j] -z[i]);
              ++node;
            }
        }
      else if (v.size () == 4u * msh.nelements)
        {
          double vloc[4];
          int node = 0;
          for (int i = 0; i < 4; ++i)
            vloc[i] = v[inode[i]];
          for (int i = 0; i < 3; ++i)
            for (int j = i + 1 ; j < 4; ++j)
              {
                vnodes[node] = vloc[j] - vloc[i];
                ++node;
              }
        }
      else
        {
          std::cerr << "bim3a_advection_upwind: "
                    << "parameter v has wrong dimension" << std::endl;
          exit (-1);
        }

      double bp[6];
      double bm[6];
      for (int i = 0; i < 6; ++i)
        {
          bp[i] = - (vnodes[i] - fabs (vnodes[i])) / 2;
          bm[i] = (vnodes[i] + fabs (vnodes[i])) / 2;
        }
      int node = 0;
      for (int ii = 0; ii < 3; ++ii)
        for (int jj = ii + 1; jj < 4; ++jj)
          {
            Lloc[ii + 4 * jj] *= bp[node];
            Lloc[jj + 4 * ii] *= bm[node];
            Lloc[5 * jj] -= Lloc[ii + 4 * jj];
            Lloc[5 * ii] -= Lloc[jj + 4 * ii];
            ++node;
          }

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          UP[inode[ii]][inode[jj]] += Lloc[ii + 4 * jj];
    }
};

void
bim3a_osc_laplacian (mesh& msh,
                     const std::vector<double>& acoeff,
                     sparse_matrix& SG)
{

  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16], p[12];
  int iel, elnodes[4];


  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          elnodes[ii] = msh.t (ii, iel);
          p[3 * ii] = msh.p (0, elnodes[ii]);
          p[3 * ii + 1] = msh.p (1, elnodes[ii]);
          p[3 * ii + 2] = msh.p (2, elnodes[ii]);
        }

      // Compute local laplacian matrix
      // and assemble into global matrix
      memset (Lloc, 0, 16 * sizeof (double));
      bim3a_osc_local_laplacian (& (msh.shg (0, 0, iel)),
                                 p, msh.volume (iel),
                                 acoeff[iel], Lloc);

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          {
            SG[elnodes[ii]][elnodes[jj]] += Lloc[ii + 4 * jj];
          }
    }

};

void
bim3a_osc_laplacian_anisotropic (mesh& msh,
                     const std::vector<double>& acoeff,
                     sparse_matrix& SG)
{

  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16], p[12];
  int iel, elnodes[4];


  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          elnodes[ii] = msh.t (ii, iel);
          p[3 * ii] = msh.p (0, elnodes[ii]);
          p[3 * ii + 1] = msh.p (1, elnodes[ii]);
          p[3 * ii + 2] = msh.p (2, elnodes[ii]);
        }

      // Compute local laplacian matrix
      // and assemble into global matrix
      memset (Lloc, 0, 16 * sizeof (double));
      bim3a_osc_local_laplacian_anisotropic (& (msh.shg (0, 0, iel)),
                                             p, msh.volume (iel),
                                             &acoeff[3 * iel], Lloc);

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          {
            SG[elnodes[ii]][elnodes[jj]] += Lloc[ii + 4 * jj];
          }
    }

};

void
bim3a_osc_advection_diffusion (mesh& msh,
                               const std::vector<double>& acoeff,
                               const std::vector<double>& v,
                               sparse_matrix& SG)
{

  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16], p[12], vloc[4];
  int iel, elnodes[4];


  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          elnodes[ii] = msh.t (ii, iel);
          p[3 * ii] = msh.p (0, elnodes[ii]);
          p[3 * ii + 1] = msh.p (1, elnodes[ii]);
          p[3 * ii + 2] = msh.p (2, elnodes[ii]);
          vloc[ii] = v[elnodes[ii]];
        }

      // Compute local laplacian matrix
      // and assemble into global matrix
      memset (Lloc, 0, 16 * sizeof (double));
      bim3a_osc_local_laplacian (& (msh.shg (0, 0, iel)),
                                 p, msh.volume (iel),
                                 acoeff[iel], Lloc);
      bim3a_local_advection (vloc, 1.0, 1.0, Lloc);

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          SG[elnodes[ii]][elnodes[jj]] += Lloc[ii + 4 * jj];

    }
};

void
bim3a_osc_advection_diffusion_anisotropic (mesh& msh,
                      const std::vector<double>& acoeff,
                      const std::vector<double>& v,
                      sparse_matrix& SG)
{
  if (SG.size () < size_t (msh.nnodes))
    SG.resize (msh.nnodes);

  double Lloc[16], p[12], vloc[4];
  int iel, elnodes[4];

  for (iel = 0; iel < msh.nelements; ++iel)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          elnodes[ii] = msh.t (ii, iel);
          p[3 * ii] = msh.p (0, elnodes[ii]);
          p[3 * ii + 1] = msh.p (1, elnodes[ii]);
          p[3 * ii + 2] = msh.p (2, elnodes[ii]);
          vloc[ii] = v[elnodes[ii]];
        }

      // Compute local laplacian matrix
      // and assemble into global matrix
      memset (Lloc, 0, 16 * sizeof (double));
      bim3a_osc_local_laplacian_anisotropic (& (msh.shg (0, 0, iel)),
                                             p, msh.volume (iel),
                                             &acoeff[3 * iel], Lloc);
      bim3a_local_advection (vloc, 1.0, 1.0, Lloc);

      for (int ii = 0; ii < 4; ++ii)
        for (int jj = 0; jj < 4; ++jj)
          SG[elnodes[ii]][elnodes[jj]] += Lloc[ii + 4 * jj];

    }
};

void
bim3a_local_advection (const double vloc[4],
                       const double a,
                       const double epsilon,
                       double lloc[16])
{
  if (epsilon != 1.0)
    for (double *ii = lloc; ii < lloc+16; ++ii)
      *(ii) *= epsilon;

  double bm, bp;
  for (double *ii = lloc; ii < lloc+16; ii+=5)
    *(ii) =  0.0;

  for (int ii = 0; ii < 3; ++ii)
    for (int jj = ii + 1; jj < 4; ++jj)
      {
        bimu_bernoulli (a * (vloc[jj] - vloc[ii]), bp, bm);
        lloc[ii + 4 * jj] *= bp;
        lloc[jj + 4 * ii] *= bm;
        lloc[5 * jj] -= lloc[ii + 4 * jj];
        lloc[5 * ii] -= lloc[jj + 4 * ii];
      }
};

void
bim3a_local_advection_jacobian (const double vloc[4],
                                const double d[4],
                                const double a,
                                const double epsilon,
                                double lloc[16])
{
  if (epsilon != 1.0)
    for (double *ii = lloc; ii < lloc+16; ++ii)
      *(ii) *= epsilon;

  double bkip, bikp, tmp;
  for (double *ii = lloc; ii < lloc+16; ii+=5)
    *(ii) =  0.0;

  for (int ii = 0; ii < 3; ++ii)
    for (int kk = ii + 1; kk < 4; ++kk)
      {
        bimu_bernoulli_derivative (a * (vloc[kk] - vloc[ii]), bkip, bikp);
        tmp =  a * (bikp * d[ii] + bkip * d[kk]);
        lloc[ii + 4 * kk] *= tmp;
        lloc[kk + 4 * ii] *= tmp;
        lloc[5 * kk] -= lloc[ii + 4 * kk];
        lloc[5 * ii] -= lloc[kk + 4 * ii];
      }
};


void
bim3a_local_laplacian (const double shg[12],
                       const double vol,
                       const double acf,
                       double Lloc[16])
{
  for (int dir = 0; dir < 3; ++dir)
    for (int ii = 0; ii < 4; ++ii)
      for (int jj = 0; jj < 4; ++jj)
    {
      Lloc[ii + 4*jj] += acf * vol *
        shg[dir + 3 * ii] * shg[dir + 3 * jj];
    }
};

void
bim3a_local_laplacian_anisotropic (const double shg[12],
                                   const double vol,
                                   const double acf[9],
                                   double Lloc[16])
{
 for (int idir = 0; idir < 3; ++idir)
  for (int jdir = 0; jdir < 3; ++jdir)
    for (int ii = 0; ii < 4; ++ii)
      for (int jj = 0; jj < 4; ++jj)
    {
      Lloc[ii + 4 * jj] += acf[idir + 3 * jdir] * vol *
        shg[idir + 3 * ii] * shg[jdir + 3 * jj];
    }
};

void
bim3a_local_laplacian_anisotropic_diag (const double shg[12],
                                        const double vol,
                                        const double acf[3],
                                        double Lloc[16])
{
 for (int dir = 0; dir < 3; ++dir)
    for (int ii = 0; ii < 4; ++ii)
      for (int jj = 0; jj < 4; ++jj)
    {
      Lloc[ii + 4*jj] += acf[dir] * vol *
        shg[dir + 3 * ii] * shg[dir + 3 * jj];
    }
};

void
bim3a_osc_local_laplacian (const double shg[12],
                           const double pts[12],
                           const double volume,
                           const double acoeff,
                           double Lloc[16])
{
  int inode, idir;
  double A[12] = {0}, Ann[4]= {0}, AidotAj[6]={0}, r[12] = {0};
  double epsilonareak  = acoeff / volume / 48.0;

  for (inode = 0; inode < 4; ++inode)
    for (idir = 0; idir < 3; ++idir)
    {
      A[idir + 3 * inode] = 3.0 * volume * shg[idir + 3 * inode];
      Ann[inode] += (A[idir + 3 * inode] * A[idir + 3 * inode]);
    }

  memset (AidotAj, 0.0, 6 * sizeof (double));
  memset (r, 0.0, 12 * sizeof (double));

  for (idir = 0; idir < 3; ++idir)
    {
      r[0]  +=                  //rik dot rjk
        (pts[idir + 3 * 2] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 2] - pts[idir + 3 * 1]);
      r[1]  +=                  //ril dot rjl
        (pts[idir + 3 * 3] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 3] - pts[idir + 3 * 1]);
      r[2]  +=                  //rij dot rkj
        (pts[idir + 3 * 1] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 1] - pts[idir + 3 * 2]);
      r[3]  +=                  //ril dot rkl
        (pts[idir + 3 * 3] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 3] - pts[idir + 3 * 2]);
      r[4]  +=                  //rij dot rlj
        (pts[idir + 3 * 1] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 1] - pts[idir + 3 * 3]);
      r[5]  +=                  //rik dot rlk
        (pts[idir + 3 * 2] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 2] - pts[idir + 3 * 3]);
      r[6]  +=                  //rji dot rki
        (pts[idir + 3 * 0] - pts[idir + 3 * 1]) *
        (pts[idir + 3 * 0] - pts[idir + 3 * 2]);
      r[7]  +=                  //rjl dot rkl
        (pts[idir + 3 * 3] - pts[idir + 3 * 1]) *
        (pts[idir + 3 * 3] - pts[idir + 3 * 2]);
      r[8]  +=                  //rji dot rli
        (pts[idir + 3 * 0] - pts[idir + 3 * 1]) *
        (pts[idir + 3 * 0] - pts[idir + 3 * 3]);
      r[9]  +=                  //rjk dot rlk
        (pts[idir + 3 * 2] - pts[idir + 3 * 1]) *
        (pts[idir + 3 * 2] - pts[idir + 3 * 3]);
      r[10] +=                  //rki dot rli
        (pts[idir + 3 * 0] - pts[idir + 3 * 2]) *
        (pts[idir + 3 * 0] - pts[idir + 3 * 3]);
      r[11] +=                  //rkj dot rlj
        (pts[idir + 3 * 1] - pts[idir + 3 * 2]) *
        (pts[idir + 3 * 1] - pts[idir + 3 * 3]);


      AidotAj[0] += A[idir + 3 * 2] * // Ak dot Al
                    A[idir + 3 * 3];

      AidotAj[1] += A[idir + 3 * 1] * // Aj dot Al
                    A[idir + 3 * 3];

      AidotAj[2] += A[idir + 3 * 1] * // Aj dot Ak
                    A[idir + 3 * 2];

      AidotAj[3] += A[idir + 3 * 0] * // Ai dot Al
                    A[idir + 3 * 3];

      AidotAj[4] += A[idir + 3 * 0] * // Ai dot Ak
                    A[idir + 3 * 2];

      AidotAj[5] += A[idir + 3 * 0] * // Ai dot Aj
                    A[idir + 3 * 1];
    }

  double tmp;

  tmp = - epsilonareak * (2.0 * r[0] * r[1] + AidotAj[0] *
      (r[0]  * r[0]  / Ann[3] + r[1] * r[1]  / Ann[2]));
  Lloc[0 + 4 * 1] += tmp;
  Lloc[1 + 4 * 0] += tmp;
  Lloc[0 + 4 * 0] -= tmp;
  Lloc[1 + 4 * 1] -= tmp;

  tmp = - epsilonareak * (2.0 * r[2] * r[3] + AidotAj[1] *
      (r[2]  * r[2]  / Ann[3] + r[3] * r[3]  / Ann[1]));
  Lloc[0 + 4 * 2] += tmp;
  Lloc[2 + 4 * 0] += tmp;
  Lloc[0 + 4 * 0] -= tmp;
  Lloc[2 + 4 * 2] -= tmp;

  tmp = - epsilonareak * (2.0 * r[4] * r[5] + AidotAj[2] *
      (r[4]  * r[4]  / Ann[2] + r[5] * r[5]  / Ann[1]));
  Lloc[0 + 4 * 3] += tmp;
  Lloc[3 + 4 * 0] += tmp;
  Lloc[0 + 4 * 0] -= tmp;
  Lloc[3 + 4 * 3] -= tmp;

  tmp = - epsilonareak * (2.0 * r[6] * r[7] + AidotAj[3] *
      (r[6]  * r[6]  / Ann[3] + r[7] * r[7]  / Ann[0]));
  Lloc[1 + 4 * 2] += tmp;
  Lloc[2 + 4 * 1] += tmp;
  Lloc[1 + 4 * 1] -= tmp;
  Lloc[2 + 4 * 2] -= tmp;


  tmp = - epsilonareak * (2.0 * r[8] * r[9] + AidotAj[4] *
      (r[8]  * r[8]  / Ann[2] + r[9] * r[9]  / Ann[0]));
  Lloc[1 + 4 * 3] += tmp;
  Lloc[3 + 4 * 1] += tmp;
  Lloc[1 + 4 * 1] -= tmp;
  Lloc[3 + 4 * 3] -= tmp;

  tmp = - epsilonareak * (2.0 * r[10] * r[11] + AidotAj[5] *
      (r[10] * r[10] / Ann[1] + r[11] * r[11] / Ann[0]));
  Lloc[2 + 4 * 3] += tmp;
  Lloc[3 + 4 * 2] += tmp;
  Lloc[2 + 4 * 2] -= tmp;
  Lloc[3 + 4 * 3] -= tmp;

};

void
bim3a_osc_local_laplacian_anisotropic (const double shg[12],
                                       const double pts[12],
                                       const double volume,
                                       const double acoeff[3],
                                       double Lloc[16])
{
  int inode, idir;
  double A[12] = {0}, Ann[4]= {0},
         AidotAj[6]={0}, r[12] = {0}, dcoeff[3]={0};
  double d = cbrt (acoeff[0] * acoeff[1] * acoeff[2]);
  double epsilonareak  = d / volume / 48.0;

  for (inode = 0; inode < 4; ++inode)
    for (idir = 0; idir < 3; ++idir)
    {
      A[idir + 3 * inode] = 3.0 * volume * shg[idir + 3 * inode];
      Ann[inode] += (A[idir + 3 * inode] * A[idir + 3 * inode]);
      dcoeff[idir] = sqrt (d / acoeff[idir]);
    }

  memset (AidotAj, 0.0, 6 * sizeof (double));
  memset (r, 0.0, 12 * sizeof (double));

  for (idir = 0; idir < 3; ++idir)
    {
      r[0]  +=                  //rik dot rjk
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 2] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 2] - pts[idir + 3 * 1]);
      r[1]  +=                  //ril dot rjl
        dcoeff[idir] * dcoeff[idir]*
        (pts[idir + 3 * 3] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 3] - pts[idir + 3 * 1]);
      r[2]  +=                  //rij dot rkj
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 1] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 1] - pts[idir + 3 * 2]);
      r[3]  +=                  //ril dot rkl
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 3] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 3] - pts[idir + 3 * 2]);
      r[4]  +=                  //rij dot rlj
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 1] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 1] - pts[idir + 3 * 3]);
      r[5]  +=                  //rik dot rlk
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 2] - pts[idir + 3 * 0]) *
        (pts[idir + 3 * 2] - pts[idir + 3 * 3]);
      r[6]  +=                  //rji dot rki
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 0] - pts[idir + 3 * 1]) *
        (pts[idir + 3 * 0] - pts[idir + 3 * 2]);
      r[7]  +=                  //rjl dot rkl
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 3] - pts[idir + 3 * 1]) *
        (pts[idir + 3 * 3] - pts[idir + 3 * 2]);
      r[8]  +=                  //rji dot rli
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 0] - pts[idir + 3 * 1]) *
        (pts[idir + 3 * 0] - pts[idir + 3 * 3]);
      r[9]  +=                  //rjk dot rlk
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 2] - pts[idir + 3 * 1]) *
        (pts[idir + 3 * 2] - pts[idir + 3 * 3]);
      r[10] +=                  //rki dot rli
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 0] - pts[idir + 3 * 2]) *
        (pts[idir + 3 * 0] - pts[idir + 3 * 3]);
      r[11] +=                  //rkj dot rlj
        dcoeff[idir] * dcoeff[idir] *
        (pts[idir + 3 * 1] - pts[idir + 3 * 2]) *
        (pts[idir + 3 * 1] - pts[idir + 3 * 3]);

      AidotAj[0] += A[idir + 3 * 2] * // Ak dot Al
                    A[idir + 3 * 3];

      AidotAj[1] += A[idir + 3 * 1] * // Aj dot Al
                    A[idir + 3 * 3];

      AidotAj[2] += A[idir + 3 * 1] * // Aj dot Ak
                    A[idir + 3 * 2];

      AidotAj[3] += A[idir + 3 * 0] * // Ai dot Al
                    A[idir + 3 * 3];

      AidotAj[4] += A[idir + 3 * 0] * // Ai dot Ak
                    A[idir + 3 * 2];

      AidotAj[5] += A[idir + 3 * 0] * // Ai dot Aj
                    A[idir + 3 * 1];
    }

  double tmp;

  tmp = - epsilonareak * (2.0 * r[0] * r[1] + AidotAj[0] *
      (r[0]  * r[0]  / Ann[3] + r[1] * r[1]  / Ann[2]));
  Lloc[0 + 4 * 1] += tmp;
  Lloc[1 + 4 * 0] += tmp;
  Lloc[0 + 4 * 0] -= tmp;
  Lloc[1 + 4 * 1] -= tmp;

  tmp = - epsilonareak * (2.0 * r[2] * r[3] + AidotAj[1] *
      (r[2]  * r[2]  / Ann[3] + r[3] * r[3]  / Ann[1]));
  Lloc[0 + 4 * 2] += tmp;
  Lloc[2 + 4 * 0] += tmp;
  Lloc[0 + 4 * 0] -= tmp;
  Lloc[2 + 4 * 2] -= tmp;

  tmp = - epsilonareak * (2.0 * r[4] * r[5] + AidotAj[2] *
      (r[4]  * r[4]  / Ann[2] + r[5] * r[5]  / Ann[1]));
  Lloc[0 + 4 * 3] += tmp;
  Lloc[3 + 4 * 0] += tmp;
  Lloc[0 + 4 * 0] -= tmp;
  Lloc[3 + 4 * 3] -= tmp;

  tmp = - epsilonareak * (2.0 * r[6] * r[7] + AidotAj[3] *
      (r[6]  * r[6]  / Ann[3] + r[7] * r[7]  / Ann[0]));
  Lloc[1 + 4 * 2] += tmp;
  Lloc[2 + 4 * 1] += tmp;
  Lloc[1 + 4 * 1] -= tmp;
  Lloc[2 + 4 * 2] -= tmp;


  tmp = - epsilonareak * (2.0 * r[8] * r[9] + AidotAj[4] *
      (r[8]  * r[8]  / Ann[2] + r[9] * r[9]  / Ann[0]));
  Lloc[1 + 4 * 3] += tmp;
  Lloc[3 + 4 * 1] += tmp;
  Lloc[1 + 4 * 1] -= tmp;
  Lloc[3 + 4 * 3] -= tmp;

  tmp = - epsilonareak * (2.0 * r[10] * r[11] + AidotAj[5] *
      (r[10] * r[10] / Ann[1] + r[11] * r[11] / Ann[0]));
  Lloc[2 + 4 * 3] += tmp;
  Lloc[3 + 4 * 2] += tmp;
  Lloc[2 + 4 * 2] -= tmp;
  Lloc[3 + 4 * 3] -= tmp;

};

void
bim3a_local_reaction (const double shp[16],
                      const double wjacdet[4],
                      const double e,
                      const double n[4],
                      double Lloc[16])
{
  for (int ii = 0; ii < 4; ++ii)
    Lloc[5*ii] += n[ii] * e * wjacdet[ii];
};

void
bim3a_local_rhs (const double shp[16],
                 const double wjacdet[4],
                 const double e,
                 const double n[4],
                 double bLoc[4])
{
  for (int ii = 0; ii < 4; ++ii)
    bLoc[ii] += n[ii] * e * wjacdet[ii];
};

void
bimu_bernoulli (double x, double &bp, double &bn)
{
  const double xlim = 1.0e-2;
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

void
bimu_bernoulli_derivative (double x, double &bpp, double &bnp)
{
  const double xlim = 1.0e-5;
  double ax  = fabs (x);

  double bp  = 0.0;
  double bn  = 0.0;
  bimu_bernoulli (x, bp, bn);

  bpp  = 0.0;
  bnp  = 0.0;

  //  X=0
  if (x == 0.0)
    {
      bpp = -.5;
      bnp =  .5;
      return;
    }

  // INTERMEDIATE VALUES
  if (ax > xlim)
    {
      bpp = (bp / x) * (1 - bn);
      bnp = - (bn / x) * (1 - bp);
      return;
    }

  // SMALL VALUES
  if (ax <= xlim &&  ax != 0.0)
    {
      bpp = -.5 + x / 6.0 - pow (x, 3) / 180.0;
      bnp =  .5 + x / 6.0 - pow (x, 3) / 180.0;
      return;
    }

};

void
bim3a_dirichlet_bc (sparse_matrix& M,
                    std::vector<double>& b,
                    const std::vector<int>& bnodes,
                    const std::vector<double>& vnodes)
{
  sparse_matrix::col_iterator j;
  for (unsigned int it = 0; it < bnodes.size (); ++it)
    {
      int i = bnodes[it];
      b[i] = vnodes[it];
      if (M[i].size ())
        for (j = M[i].begin (); j != M[i].end (); ++j)
          {
            int jj = M.col_idx (j);
            if (jj != i)
              {
                M[i][jj] = 0.0;
                b[jj] -= M[jj][i] * b[i];
                M[jj][i] = 0.0;
              }
          }
      b[i] *= M[i][i];
    }
}

void
bim3a_pde_gradient (mesh& msh,
                   const std::vector<double>& u,
                   std::vector<double>& g)
{
  g.clear ();
  g.resize (msh.nelements * 3);
  for (int iel = 0; iel < msh.nelements; ++iel)
    for (int inode = 0; inode < 4; ++inode)
      {
        g[iel * 3 + 0] +=
          msh.shg (0, inode, iel) * u[msh.t (inode, iel)];
        g[iel * 3 + 1] +=
          msh.shg (1, inode, iel) * u[msh.t (inode, iel)];
        g[iel * 3 + 2] +=
          msh.shg ( 2, inode, iel) * u[msh.t (inode, iel)];
      }
}

void
bim3a_norm (mesh& msh,
           const std::vector<double>& v,
           double& norm,
           norm_type type)
{
  if (type == Inf)
    {
      norm = 0.0;
      for (unsigned int i = 0; i < v.size (); ++i)
        {
          double temp = fabs (v[i]);
          if (norm < temp)
            norm = temp;
        }
    }
  else if (type == L2 || type == H1)
    {
      norm = 0.0;
      std::vector<double> ecoeff (msh.nelements, 1.0);
      std::vector<double> ncoeff (msh.nnodes, 1.0);
      sparse_matrix M;
      bim3a_reaction (msh, ecoeff, ncoeff, M);
      if (type == H1)
        bim3a_laplacian (msh, ecoeff, M);
      std::vector<double> temp;
      temp = M * v;
      for (unsigned int i = 0; i < v.size (); ++i)
        norm += v[i] * temp[i];
      norm = sqrt (norm);
    }
}

//}
