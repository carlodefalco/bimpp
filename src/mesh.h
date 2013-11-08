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


/*! \file mesh.h
  \brief Classes and methods for managing tetrahedral meshes.
*/

#ifndef HAVE_MESH_H
#define HAVE_MESH_H 1

#include <string>

//namespace bim
//{
/// Tetrahedral FEM grid.
class mesh 
{

protected:

  /// Points array.
  double *p_data;
  /// Elements (tetrahedra)  array.
  int    *t_data;
  /// Boundary faces array.
  int    *e_data;

  double *shp_data,      //!< Shape function values array.
         *shg_data,      //!< Shape function gradients values array.
         *wjacdet_data,  //!< weighted Jacobian det values array.
         *volume_data;   //!< Tetrahedra volume values array.

  mesh () { };
  
public:
  
  int nnodes,    //!< Number of mesh nodes.
      nelements, //!< Number of mesh tetrahedra.
      nfaces;    //!< Number of (boundary?) faces

  void read (std::string filename); //!< Reads a mesh file.
  void write (std::string filename); //!< Writes a mesh file.

  mesh (std::string filename) :  p_data (NULL),
    t_data (NULL),
    e_data (NULL),
    shp_data (NULL),
    shg_data (NULL),
    wjacdet_data (NULL),
    volume_data (NULL)
      {read (filename);} //!< Build a mesh by reading from file.

  virtual ~mesh ();     //!< Destructor (frees the pointer members).

  void precompute_properties (); //!< Cache some mesh properties 
                                 //! required by the discrete operator constructors
                                 //! (shape functions, volumes, etc.).

  /// i-th coordinate of a node.
  double inline 
    &p (int idir, int inode) 
    {return (*(p_data+idir+3*inode));}; 

  /// i-th node of a tetrahedra.
  int inline 
    &t (int inode, int iel) 
    {return (*(t_data+inode+5*iel));};

  /// i-th boundary face
  int inline 
    &e (int ient, int ifc) 
    {return (*(e_data+ient+10*ifc));};

  /// i-th coordinate of a node (const version).
  const double inline 
    &p (int idir, int inode) const 
    {return (*(p_data+idir+3*inode));};
  
  /// i-th node of a tetrahedra (const version).
  const int inline 
    &t (int inode, int iel) const 
    {return (*(t_data+inode+5*iel));};
  
  /// i-th boundary face (const version).
  const int inline 
    &e (int ient, int ifc) const 
    {return (*(e_data+ient+10*ifc));};

  /// node value of a node shape function.
  double inline 
    &shp (int inode, int jnode) 
    {return (*(shp_data+inode+4*jnode));};
  
  /// value of node shape function gradient component on an element.
  double inline 
    &shg (int idir, int inode, int iel) 
    {return (*(shg_data+idir+3*(inode+(4*iel))));};
  
  /// determinant of the jacobian of a shape function on an element.
  double inline 
    &wjacdet (int inode, int iel) 
    {return (*(wjacdet_data+inode+4*iel));};
  
  /// volume of an element.
  double inline 
    &volume (int iel) 
    {return (*(volume_data+iel));};


  /// node value of a node shape function.
  const double inline 
    &shp (int inode, int jnode) const 
    {return (*(shp_data+inode+4*jnode));};
  
  /// value of node shape function gradient component on an element.
  const double inline 
    &shg (int idir, int inode, int iel) const 
    {return (*(shg_data+idir+3*(inode+(4*iel))));};
  
  /// determinant of the jacobian of a shape function on an element.
  const double inline 
    &wjacdet (int inode, int iel) const 
    {return (*(wjacdet_data+inode+4*iel));};
  
  /// volume of an element.
  const double inline 
    &volume (int iel) const 
    {return (*(volume_data+iel));};

  /// mesh textual output.
  friend std::ostream &operator<< (std::ostream &, mesh &);
};

//}

#endif

