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

  using idx_t = p4est_gloidx_t;

  // forward declaration of friend class
  class  quadrant_t;

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
    
    /// Operator -> to get access to the wrapped pointer.
    quadrant_t *
    operator-> ()
    { return this->data; };
    
    /// Operator ->. const version.
    const quadrant_t *
    operator-> () const
    { return this->data; };

    /// Compare two quadrant_iterator objects.
    bool
    operator== (const quadrant_iterator& other)
    { return (this->data == other.data); };

    /// Return the opposite of ==.
    bool
    operator!= (const quadrant_iterator& other)
    { return !((*this) == other); };

    /// Default constructor.
    quadrant_iterator (quadrant_t *_data = nullptr) :
      data (_data)
    { };

    /// Move to first forest quadrant
    void
    reset ();

  private:
    quadrant_t *data;      
  };

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
    { };

    /// Get the i-th coordinate of the j-th vertex.
    double
    p (idx_t i, idx_t j);

    /// Get rank-local index of the i-th vertex
    idx_t
    t (idx_t i);

    /// Get global index of the i-th vertex
    idx_t
    gt (idx_t i);

    /// True if the i-th vertex is hanging.
    bool
    is_hanging (idx_t i);

    /// Return the ip-th parent for the in-th vertex.
    int
    parent (idx_t ip, idx_t in);
      
    /// Index of the edge of the current tree
    //  on which the i-th vertex lies, return
    //  NOT_ON_BOUNDARY if an interior vertex.
    static const idx_t NOT_ON_BOUNDARY = P4EST_ROOT_LEN + 1;
    idx_t
    e (idx_t i);
    
    /// Return index of current quadrant across all trees on current process.
    p4est_locidx_t
    get_forest_quad_idx ()
    { return forest_quad_idx; };

    /// Return index of current quadrant across all trees on current process.
    p4est_gloidx_t
    get_global_quad_idx (MPI_Comm comm = MPI_COMM_WORLD)
    {
      int rank;
      MPI_Comm_rank (comm, &rank);
      return forest_quad_idx +
        the_tmesh->p4est->global_first_quadrant[rank];
    };

    /// Return index of current tree.
    p4est_locidx_t
    get_tree_idx ()
    { return tree_idx; };
    
    /// Update stored data.
    void
    update (p4est_topidx_t tree,
            p4est_quadrant_t *q);
    
    /// A pointer to the owning mesh is needed to get physical mapping. 
    tmesh               *the_tmesh; 
    p4est_tree_t        *tree;
    p4est_quadrant_t    *the_quadrant;

    friend class tmesh::quadrant_iterator;
    
  private:

    sc_array_t           *tquadrants;
    p4est_locidx_t        num_quadrants;    // Q
    // Local and global indices for looping.
    p4est_topidx_t        tree_idx;         // tt
    p4est_locidx_t        forest_quad_idx;  // k 
    p4est_locidx_t        tree_quad_idx;    // q

    /// Buffer used when quering coordinates.
    double                vxyz[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    idx_t                 tbuff[4] = {0,0,0,0};
    bool                  hbuff[4] = {false,false,false,false};
    int                   pbuff[8] = {-1,-1,-1,-1,-1,-1,-1,-1};
  };

  /// Default constructor, set all pointers to nullptr.
  tmesh ()
    : p4est (nullptr), conn (nullptr),
      current_quadrant (this, 0, nullptr),
      lnodes (nullptr), mesh (nullptr)
  { };

  /// Load a p4est and connectivity from a file.
  tmesh (const char *filename)
    : p4est (nullptr), current_quadrant (this, 0, nullptr),
      lnodes (nullptr), mesh (nullptr)
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

  /// Export nodal field f to a octbin.gz file for visualization.
  void
  octbin_export (const char * filename,
                 const std::vector<double> & f,
                 MPI_Comm comm = MPI_COMM_WORLD);
  
  /// Get an iterator to the first quadrant of the mesh.
  quadrant_iterator
  begin_quadrant_sweep ();

  /// Get a null quadrant iterator to signal end of the sweep.
  quadrant_iterator
  end_quadrant_sweep ()
  { return quadrant_iterator (); };

  /// Get an iterator to the first neighbor of the current quadrant.
  quadrant_iterator
  begin_neighbor_sweep ();

  /// Get a null quadrant iterator to signal end of the sweep.
  quadrant_iterator
  end_neighbor_sweep ()
  { return quadrant_iterator (); };

  /// Set functor to mark quadrants for refinement.
  void
  set_refine_marker
  (std::function<int (quadrant_iterator)> fun)
  { refine_marker = fun; };

  /// Set functor to mark quadrants for coarsening.
  void
  set_coarsen_marker
  (std::function<int (std::function<void (idx_t, quadrant_iterator&)>)> fun)
  { coarsen_marker = fun; };

  /// Refine marked quadrants, balance the quadtree and
  /// re-partition over the processors.
  void
  refine (int recursive = 0, int partforcoarsen = 1);

  /// Coarsen marked quadrants, balance the quadtree and
  /// re-partition over the processors.
  void
  coarsen (int recursive = 0, int partforcoarsen = 1);

  /// Compute nodes numbering and quadrant neighbours.
  void
  update ();

  /// Return number of nodes owned by local process
  idx_t
  num_owned_nodes ()
  {
    if (! lnodes) update ();
    return lnodes->owned_count;
  };

  /// Return number of nodes of quadrants owned or shared by local process
  idx_t
  num_local_nodes ()    
  {
      if (! lnodes) update ();
      return lnodes->num_local_nodes;
  };

  /// Return total number of quadrants owned by all process
  idx_t
  num_global_nodes (MPI_Comm comm)
  {
    idx_t retval = 0;
    int size;
    MPI_Comm_size (comm, &size);
    if (! lnodes) update ();
    for (int i = 0; i < size; ++i)
      retval += lnodes->global_owned_count[i];
    return retval;
  };

  /// Return number of quadrants owned by local process
  idx_t
  num_local_quadrants ()    
  {
      if (! lnodes) update ();
      return lnodes->num_local_elements;
  };
  
  /// P4EST pointers describing the tmesh,
  /// temporarily public untli the API is stable.
  p4est_t              *p4est;
  p4est_connectivity_t *conn;
  quadrant_t            current_quadrant;
  p4est_lnodes_t       *lnodes;
  p4est_mesh_t         *mesh;
  
private:
  
  std::function<int (quadrant_iterator)> refine_marker;
  std::function<int (std::function<void (idx_t, quadrant_iterator&)>)> coarsen_marker;

  static int
  refine_callback (p4est_t*, p4est_topidx_t, p4est_quadrant_t*);

  static int
  coarsen_callback (p4est_t*, p4est_topidx_t, p4est_quadrant_t* []);
  
  static void
  select_quad (tmesh *_tmesh,
               p4est_topidx_t tree_idx,
               p4est_quadrant_t* qt [],
               idx_t ii, quadrant_iterator& qi);
  
};


#endif /* TMESH_H */

