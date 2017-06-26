/*
  Copyright (C) 2017 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file tmesh.h
  \brief Interface for p4est library
*/

#ifndef TMESH_H
#define TMESH_H


#include <octave_file_io.h>

#include <mpi.h>

#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>

#include <cassert>
#include <functional>
#include <array>
#include <vector>

/// C++ interface class for p4est 2d quadrant meshes.
class
tmesh
{

public:

  using idx_t = unsigned int;

  /// C++ interface class to access properties of the
  /// current quadrant.
  class
  quadrant_t
  {

  public:

    /// Simple constructor needs at least a pointer
    /// to the container tmesh.
    quadrant_t (tmesh *_tmesh,
                p4est_topidx_t _tree = 0,
                p4est_quadrant_t *_quadrant = nullptr) :
      the_tmesh(_tmesh), tree_idx(_tree), the_quadrant(_quadrant)
    {  };

    /// Get the i-th coordinate of the j-th vertex.
    double
    p (idx_t i, idx_t j);

    /// Get global index of the i-th vertex
    idx_t
    t (idx_t i);

    /// True if the i-th vertex is hanging.
    bool
    is_hanging (idx_t i);

    /// Return the list of parents for a hanging vertex.
    void
    get_parents (idx_t i, std::vector<idx_t> &pv);

    /// Index of the boundary side on which the i-th vertex lies,
    /// 0 for interior vertices.
    idx_t
    e (idx_t i);

    /// Update stored data.
    void
    update (p4est_topidx_t tree,
            p4est_quadrant_t *q);

    /// A pointer to the owning mesh is needed to get physical mapping. 
    tmesh               *the_tmesh; 

    p4est_tree_t         *tree;
    sc_array_t           *tquadrants;
    p4est_locidx_t        num_quadrants;    // Q
    // Local and global indices for looping.
    p4est_topidx_t        tree_idx;         // tt
    p4est_locidx_t        forest_quad_idx;  // k 
    p4est_locidx_t        tree_quad_idx;    // q
    p4est_quadrant_t     *the_quadrant;
    
  private:


    /// Buffer used when quering coordinates.
    double                vxyz[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
  };


  /// Iterator to sweep through the quadrants of a tmesh.
  /// This is essentially a decorator of quadrant_t*.
  class
  quadrant_iterator  
  {

  public:

    /// Get next quadrant.
    void 
    operator++ ();

    /// Dereference.
    quadrant_t&
    operator* ()
    { return *(this->data); };

    /// Dereference. const version.
    const quadrant_t&
    operator* () const
    { return *(this->data); };

    /// Get direct access to the wrapped pointer. 
    quadrant_t*
    get_data ()
    { return this->data; };

    /// Get direct access to the wrapped pointer. const version.
    const quadrant_t*
    get_data () const
    { return this->data; };

    /// Compare two quadrant_iterator objects.
    bool
    operator== (const quadrant_iterator& other)
    { return (this->get_data () == other.get_data ()); };

    /// Return the opposite of ==.
    bool
    operator!= (const quadrant_iterator& other)
    { return !((*this) == other); };

    /// Default constructor.
    quadrant_iterator (quadrant_t *_data = nullptr) :
      data (_data)
    { };

  private:
    quadrant_t *data;      
  };

  /// Default constructor, set all pointers to nullptr.
  tmesh ()
    : p4est (nullptr), conn (nullptr), current_quadrant (this, 0, nullptr)
  {  };

  /// Load a p4est and connectivity from a file.
  tmesh (const char *filename)
    : p4est (nullptr), current_quadrant (this, 0, nullptr)
  { load (filename); };


  ~tmesh ();

  /// Load a connectivity from a compressed binary
  /// Octave file then init the p4est.
  void
  read_connectivity (const char *filename,
                     int source = 0,
                     MPI_Comm comm = MPI_COMM_WORLD);

  /// Save the p4est and connectivity to a file.
  void
  save (const char *filename);

  /// Load the p4est and connectivity from a file.
  void
  load (const char *filename,
        MPI_Comm comm = MPI_COMM_WORLD);

  /// Export exploded mesh to a vtk file for visualization.
  void
  vtk_export (const char *filename);

  /// Get an iterator to the first quadrant of the mesh.
  quadrant_iterator
  begin_quadrant_sweep ();

  /// Get a null quadrant iterator to signal end of the sweep.
  quadrant_iterator
  end_quadrant_sweep ()
  { return quadrant_iterator (); };

  /// Set functor to mark quadrants for refinement.
  void
  set_refine_marker
  (std::function<int (quadrant_iterator)> fun)
  { refine_marker = fun; };

  /// Set functor to mark quadrants for coarsening.
  void
  set_coarsen_marker
  (std::function<int (quadrant_iterator[4])> fun)
  { coarsen_marker = fun; };

  /// Refine marked quadrants, balance the quadtree and
  /// re-partition over the processors.
  void
  refine (int recursive = 0, int partforcoarsen = 0);

  /// P4EST pointers describing the tmesh,
  /// temporarily public untli the API is stable.
  p4est_t              *p4est;
  p4est_connectivity_t *conn;
  quadrant_t            current_quadrant;

private:

  std::function<int (quadrant_iterator)> refine_marker;
  std::function<int (quadrant_iterator[4])> coarsen_marker;

  static int
  refine_callback (p4est_t*, p4est_topidx_t, p4est_quadrant_t*);

  static int
  coarsen_callback (p4est_t*, p4est_topidx_t, p4est_quadrant_t* []);

};


#endif /* TMESH_H */

