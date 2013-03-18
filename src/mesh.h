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

#ifndef HAVE_MESH_H
#define HAVE_MESH_H 1

#include <string>

/// Tetrahedral FEM grid.
class mesh 
{

  double *p_data;
  int    *t_data;
  int    *e_data;

  double *shp_data, *shg_data, *wjacdet_data, *volume_data;

 public:
  
  int nnodes, nelements, nfaces;

  void read (std::string filename);
  void write (std::string filename);

  mesh (std::string filename) {read (filename);}

  void precompute_properties ();

  double &p (int idir, int inode) {return (*(p_data+idir+3*inode));};
  int    &t (int inode, int iel) {return (*(t_data+inode+5*iel));};
  int    &e (int ient, int ifc) {return (*(e_data+ient+10*ifc));};

  const double &p (int idir, int inode) const {return (*(p_data+idir+3*inode));};
  const int    &t (int inode, int iel) const {return (*(t_data+inode+5*iel));};
  const int    &e (int ient, int ifc) const {return (*(e_data+ient+10*ifc));};

  double &shp (int inode, int jnode) {return (*(shp_data+inode+4*jnode));};
  double &shg (int idir, int inode, int iel) {return (*(shg_data+idir+3*(inode+(4*iel))));};
  double &wjacdet (int inode, int iel) {return (*(wjacdet_data+inode+4*iel));};
  double &volume (int iel) {return (*(volume_data+iel));};

  const double &shp (int inode, int jnode) const {return (*(shp_data+inode+4*jnode));};
  const double &shg (int idir, int inode, int iel) const {return (*(shg_data+idir+3*(inode+(4*iel))));};
  const double &wjacdet (int inode, int iel) const {return (*(wjacdet_data+inode+4*iel));};
  const double &volume (int iel) const {return (*(volume_data+iel));};

  friend std::ostream &operator<< (std::ostream &, mesh &);
};

#endif

