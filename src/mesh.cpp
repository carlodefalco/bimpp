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

/*! \file mesh.cpp
  \brief Classes and methods for managing tetrahedral meshes.
*/

#include "mesh.h"
#include <iostream>
#include <fstream>
#include <new>
#include <cstring>
#include <algorithm>
#include <cmath>

//namespace bim
//{

mesh::~mesh ()
{
  delete [] p_data;
  delete [] t_data;
  delete [] e_data;
  delete [] shp_data;
  delete [] shg_data;
  delete [] wjacdet_data;
  delete [] volume_data;
};

int mesh::read (std::string filename)
{
  int retval = 0;
  std::ifstream fin (filename.c_str ());
  if (! fin)
    {
      std::cout << "error opening file "
                << filename << std::endl;
      retval = -1;
    }
  else
    {
      fin >> nnodes >> nelements >> nfaces;
      p_data = new (std::nothrow) double[3*nnodes];
      t_data = new (std::nothrow) int[5*nelements];
      e_data = new (std::nothrow) int[10*nfaces];

      for (int ii = 0; ii < nnodes; ++ii)
        for (int kk = 0; kk < 3; ++kk)
          fin >> p (kk, ii);

      for (int ii = 0; ii < nelements; ++ii)
        for (int kk = 0; kk < 5; ++kk)
          fin >> t (kk, ii);

      for (int ii = 0; ii < nfaces; ++ii)
        for (int kk = 0; kk < 10; ++kk)
          fin >> e (kk, ii);

      fin.close ();
      retval = 0;
    }
  return retval;
};


void mesh::write (std::string filename)
{
  std::ofstream fout (filename.c_str ());
  if (! fout) std::cout << "error opening file " 
                        << filename << std::endl;
  fout << (*this);
  fout.close ();
};

void mesh::precompute_properties ()
{

  double weight[4]    = {.25, .25, .25, .25};
  double Nb2, Nb3, Nb4, detJ, Kkvolume;

  shp_data     = new (std::nothrow) double[4*4];
  if (! shp_data) std::cout << "error allocating memory" 
                            << std::endl;

  shg_data     = new (std::nothrow) double[3*4*nelements];
  if (! shg_data) std::cout << "error allocating memory" 
                            << std::endl;

  wjacdet_data = new (std::nothrow) double[4*nelements];
  if (! wjacdet_data) std::cout << "error allocating memory" 
                                << std::endl;

  volume_data  = new (std::nothrow) double[nelements];
  if (! volume_data) std::cout << "error allocating memory" 
                               << std::endl;

  //std::cout << "compute shp" << std::endl;
  for (int inode = 0; inode < 4; ++inode)
    for (int jnode = 0; jnode < 4; ++jnode)
      if (inode == jnode)
        shp (inode, jnode) = 1.0;
      else
        shp (inode, jnode) = 0.0;

  //std::cout << "start element loop" << std::endl;

  for (int iel = 0; iel < nelements; ++iel)
    {
      //std::cout << "element number " << iel << std::endl;
      double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;

      x1 = p (0, t (0, iel));
      y1 = p (1, t (0, iel));
      z1 = p (2, t (0, iel));
      x2 = p (0, t (1, iel));
      y2 = p (1, t (1, iel));
      z2 = p (2, t (1, iel));
      x3 = p (0, t (2, iel));
      y3 = p (1, t (2, iel));
      z3 = p (2, t (2, iel));
      x4 = p (0, t (3, iel));
      y4 = p (1, t (3, iel));
      z4 = p (2, t (3, iel));

      Nb2 = y1 * (z3-z4) + y3 * (z4-z1) + y4 * (z1-z3);
      Nb3 = y1 * (z4-z2) + y2 * (z1-z4) + y4 * (z2-z1);
      Nb4 = y1 * (z2-z3) + y2 * (z3-z1) + y3 * (z1-z2);

      // Determinant of the Jacobian of the
      // transformation from the base tetrahedron
      // to the tetrahedron K
      detJ = (x2-x1) * Nb2 + (x3-x1) * Nb3 + (x4-x1) * Nb4;
      //std::cout << "\tcomputed detJ = " << detJ << std::endl;

      // Volume of the reference tetrahedron
      Kkvolume = 1.0/6.0;

      for (int inode = 0; inode < 4; ++inode)
        wjacdet (inode, iel) = Kkvolume * weight[inode] * detJ;

      volume (iel) = Kkvolume * detJ;

      //std::cout << "\tcomputed wjacdet = " << wjacdet (0, iel) << std::endl;

      // Shape function gradients follow
      // first index represents space direction
      // second index represents the shape function
      // third index represents the tetrahedron number
      shg (0, 0, iel) = (y2 * (z4-z3) + y3 * (z2-z4) 
                         + y4 * (z3-z2)) / detJ;
      shg (1, 0, iel) = (x2 * (z3-z4) + x3 * (z4-z2) 
                         + x4 * (z2-z3)) / detJ;
      shg (2, 0, iel) = (x2 * (y4-y3) + x3 * (y2-y4) 
                         + x4 * (y3-y2)) / detJ;

      shg (0, 1, iel) = Nb2 / detJ;
      shg (1, 1, iel) = (x1 * (z4-z3) + x3 * (z1-z4) 
                         + x4 * (z3-z1)) / detJ;
      shg (2, 1, iel) = (x1 * (y3-y4) + x3 * (y4-y1) 
                         + x4 * (y1-y3)) / detJ;

      shg (0, 2, iel) = Nb3 / detJ;
      shg (1, 2, iel) = (x1 * (z2-z4) + x2 * (z4-z1) 
                         + x4 * (z1-z2)) / detJ;
      shg (2, 2, iel) = (x1 * (y4-y2) + x2 * (y1-y4) 
                         + x4 * (y2-y1)) / detJ;

      shg (0, 3, iel) = Nb4 / detJ;
      shg (1, 3, iel) = (x1 * (z3-z2) + x2 * (z1-z3) 
                         + x3 * (z2-z1)) / detJ;
      shg (2, 3, iel) = (x1 * (y2-y3) + x2 * (y3-y1) 
                         + x3 * (y1-y2)) / detJ;

      //std::cout << "\tcomputed shg" << std::endl;
    }
};

void
mesh::boundary_nodes (const std::vector<int>& sidelist,
                      std::vector<int>& bnodes)
{
  bnodes.clear ();
  for (int i = 0; i < nfaces; ++i)
    if (find (sidelist.begin (), sidelist.end (), e (9, i))
        != sidelist.end ())
      {
        for (int j = 0; j < 3; ++j)
          if (find (bnodes.begin (), bnodes.end (), e (j, i))
             == bnodes.end ())
           bnodes.push_back (e (j, i));
      }
  sort (bnodes.begin (), bnodes.end ());
}

std::ostream & operator<< (std::ostream &stream, mesh &msh)
{
  stream << "nnodes = " << msh.nnodes
         << "; nelements = " << msh.nelements
         << "; nfaces = " << msh.nfaces << ";" << std::endl;

  stream << "p = [ ";
  for (int inode = 0; inode < msh.nnodes; ++inode)
    {
      for (int idir = 0; idir < 3; ++idir)
        stream << msh.p (idir, inode) << " ";
      stream << ";" << std::endl;
    }
  stream << "]';" << std::endl;

  stream << "t = [ ";
  for (int iel = 0; iel < msh.nelements; ++iel)
    {
      for (int inode = 0; inode < 5; ++inode)
        stream << msh.t (inode, iel) << " ";
      stream << ";" << std::endl;
    }
  stream << "]'; t (1:4, :) += 1;" << std::endl;

  stream << "e = [ ";
  for (int ifc = 0; ifc < msh.nfaces; ++ifc)
    {
      for (int ient = 0; ient < 10; ++ient)
        stream << msh.e (ient, ifc) << " ";
      stream << ";" << std::endl;
    }
  stream << "]'; e (1:3, :) += 1;" << std::endl;

  return stream;
}

void
mesh::element_data::update (int iel_)
{
  memcpy (&iel, &iel_, sizeof (int));

  memcpy (t, m.t_data + 5 * iel, 5 * sizeof (int));

  for (int inode = 0; inode < 4; ++inode)
    memcpy (p + 3 * inode, m.p_data + 3 * t[inode], 3 * sizeof (double));

  memcpy (shg, m.shg_data + 12 * iel, 12 * sizeof (double));

  memcpy (wjacdet, m.wjacdet_data + 4 * iel, 4 * sizeof (double));

  memcpy (&volume, m.volume_data + iel, sizeof (double));
};

//}
