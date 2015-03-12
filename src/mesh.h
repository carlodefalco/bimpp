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
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>

//namespace bim
//{
/// Tetrahedral FEM grid.
class mesh
{

protected:

  /// Points.
  double *p_data;

  /// Elements (tetrahedra).
  int    *t_data;

  /// Boundary faces.
  int    *e_data;


  double *shp_data, //!< Shape function values.
    *shg_data,      //!< Shape function gradients.
    *wjacdet_data,  //!< Weighted determinant of map Jacobian.
    *volume_data;   //!< Volume of tetrahedra.



public:

  int nnodes,    //!< Number of mesh nodes.
    nelements,   //!< Number of mesh tetrahedra.
    nfaces;      //!< Number of (boundary?) faces

  int read (std::string filename); //!< Reads a mesh file.
  void write (std::string filename); //!< Writes a mesh file.

  mesh () :
    p_data (NULL),
    t_data (NULL),
    e_data (NULL),
    shp_data (NULL),
    shg_data (NULL),
    wjacdet_data (NULL),
    volume_data (NULL) { };

  mesh (std::string filename) :
    p_data (NULL),
    t_data (NULL),
    e_data (NULL),
    shp_data (NULL),
    shg_data (NULL),
    wjacdet_data (NULL),
    volume_data (NULL)
  {
    if (! read (filename) == 0) //!< Build a mesh by reading from file.
      {
        std::cerr << "could not initialize mesh" << std::endl;
        exit (-1);
      }
  }

  virtual ~mesh ();     //!< Destructor (frees the pointer members)

  void precompute_properties (); //!< Cache some mesh properties
                                 //! required by the discrete
                                 //! operator constructors
                                 //! (shape functions,
                                 //! volumes, etc.)

  /// i-th coordinate of a node.
  inline double&
  p (int idir, int inode)
  {return (*(p_data + idir + 3 * inode));};

  /// i-th node of a tetrahedron.
  inline int&
  t (int inode, int iel)
  {return (*(t_data + inode + 5 * iel));};

  /// i-th boundary face
  inline int&
  e (int ient, int ifc)
  {return (*(e_data + ient + 10 * ifc));};

  /// i-th coordinate of a node (const version).
  inline const double&
  p (int idir, int inode) const
  {return (*(p_data + idir + 3 * inode));};

  /// i-th node of a tetrahedra (const version).
  inline const int&
  t (int inode, int iel) const
  {return (*(t_data + inode + 5 * iel));};

  /// i-th boundary face (const version).
  inline const int&
  e (int ient, int ifc) const
  {return (*(e_data + ient + 10 * ifc));};

  /// node value of a node shape function.
  inline double&
  shp (int inode, int jnode)
  {return (*(shp_data + inode + 4 * jnode));};

  /// value of node shape function gradient component on an element
  inline double&
  shg (int idir, int inode, int iel)
  {return (*(shg_data + idir + 3 * (inode + (4 * iel))));};

  /// determinant of the jacobian of a shape function on an element
  inline double&
  wjacdet (int inode, int iel)
  {return (*(wjacdet_data + inode + 4 * iel));};

  /// volume of an element.
  inline double&
  volume (int iel)
  {return (*(volume_data + iel));};

  /// node value of a node shape function.
  inline const double&
  shp (int inode, int jnode) const
  {return (*(shp_data + inode + 4 * jnode));};

  /// value of node shape function gradient component on an element
  inline const double&
  shg (int idir, int inode, int iel) const
  {return (*(shg_data + idir + 3 * (inode + (4 * iel))));};

  /// determinant of the jacobian of a shape function on an element
  inline const double&
  wjacdet (int inode, int iel) const
  {return (*(wjacdet_data + inode + 4 * iel));};

  /// volume of an element.
  inline const double&
  volume (int iel) const
  {return (*(volume_data + iel));};

  /// Find nodes in sidelist boundaries.
  void
  boundary_nodes (const std::vector<int>& sidelist,
                  std::vector<int>& bnodes);

  /// mesh textual output.
  friend std::ostream &operator<< (std::ostream &, mesh &);

  //!< cache for data of a single element.
  friend class element_data;

  class element_data
  {

  public:

    union
    {
      double dbuffer[45];
      struct
      {
        double p[12];
        double shg[12];
        double shp[16];
        double wjacdet[4];
        double volume;
      };
    };

    union
    {
      int ibuffer[6];
      struct
      {
        int iel;
        int t[5];
      };
    };


    const mesh& m;
    element_data (const mesh& _m)
      : m (_m)
    {
      memcpy (shp, m.shp_data, 16 * sizeof (double));
    };

    void
    update (int iel_);

  };

};


//}

#endif
