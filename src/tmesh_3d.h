/*! \file tmesh_3d.h
  \brief Interface for p4est library
*/

#ifndef TMESH_3D_H
#define TMESH_3D_H


#include <octave_file_io.h>

#include <mpi.h>

#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_vtk.h>

#include <cassert>
#include <functional>
#include <array>
#include <vector>


/// C++ interface class for p4est 3d octant meshes.
class
tmesh_3d
{

public:

  using idx_t = p4est_gloidx_t;

  // forward declaration of friend class
  class  octant_t;

  /// Iterator to sweep through the octants of a tmesh.
  /// This is essentially a decorator of octant_t*.
  class
  octant_iterator  
  {

  public:

    /// Get next octant.
    void 
    operator++ ();
    
    /// Dereference.
    octant_t&
    operator* ()
    { return *(this->data); };
    
    /// Dereference. const version.
    const octant_t&
    operator* () const
    { return *(this->data); };
    
    /// Operator -> to get access to the wrapped pointer.
    octant_t *
    operator-> ()
    { return this->data; };
    
    /// Operator ->. const version.
    const octant_t *
    operator-> () const
    { return this->data; };

    /// Compare two octant_iterator objects.
    bool
    operator== (const octant_iterator& other)
    { return (this->data == other.data); };

    /// Return the opposite of ==.
    bool
    operator!= (const octant_iterator& other)
    { return !((*this) == other); };

    /// Default constructor.
    octant_iterator (octant_t *_data = nullptr) :
      data (_data)
    { };

    /// Move to first forest quadrant
    void
    reset ();

  protected:
    octant_t *data;      
  };
  
  /// Iterator to sweep through the octants of a tmesh.
  /// This is essentially a decorator of octant_t*.
  class
  neighbor_iterator : public octant_iterator
  {
  public:
    /// Get next neighbor.
    void
    operator++ ();
    
    /// Default constructor.
    neighbor_iterator (octant_t *_data = nullptr, int _face_idx = -1) :
      octant_iterator (_data),
      face_neighbor (new p8est_mesh_face_neighbor_t),
      face_idx (_face_idx)
    { };
    
    /// Get the face index associated to the current neighbor.
    int get_face_idx ()
    { return face_idx; };
    
    friend class tmesh_3d::octant_t;
    
  private:
    p8est_mesh_face_neighbor_t * face_neighbor; // mfn
    
    int face_idx; /// Face index in 0...5 (-1 if not defined). CESARE(controlla 0...5)
  };
  
  /// C++ interface class to access properties of the
  /// current octant.
  class
  octant_t
  {

  public:

    /// Simple constructor needs at least a pointer
    /// to the container tmesh.
    octant_t (tmesh_3d *_tmesh,
                p4est_topidx_t _tree = 0,
                p8est_quadrant_t *_octant = nullptr) :
      the_tmesh(_tmesh), the_octant(_octant), tree_idx(_tree),
      is_ghost(false), qtq(-1)
    { };

    /// Get the i-th coordinate of the j-th vertex.
    double
    p (idx_t i, idx_t j);
    
    /// Get the i-th coordinate of the centroid.
    double
    centroid (idx_t i);
    
    /// Get the i-th coordinate of the centroid of the j-th face
    double
    face_centroid (idx_t i, int j); //CESARE(cambia tipo j da int a indice di faccia globale)

    /// Get rank-local index of the i-th vertex, or global index for ghosts.
    idx_t
    t (idx_t i);

    /// Get global index of the i-th vertex
    idx_t
    gt (idx_t i);

    /// True if the i-th vertex is hanging.
    bool
    is_hanging (idx_t i);

    /// Return the rank-local (or global for ghosts)
    /// ip-th parent for the in-th vertex.
    /// TODO : type should be p4est_gloidx_t
    int
    parent (idx_t ip, idx_t in);
      
    /// Return the ip-th parent for the in-th vertex,
    /// use global numbering.
    /// TODO : type should be p4est_gloidx_t
    int
    gparent (idx_t ip, idx_t in);

    static const idx_t NOT_ON_BOUNDARY = P8EST_ROOT_LEN + 1;
    /// Index of the edge of the current tree
    //  on which the i-th vertex lies, return
    //  NOT_ON_BOUNDARY if an interior vertex.
    idx_t
    e (idx_t i);
    
    /// Get an iterator to the first neighbor
    /// of the current octant.
    neighbor_iterator
    begin_neighbor_sweep ();
    
    /// Get a null octant iterator to signal end of the sweep.
    neighbor_iterator
    end_neighbor_sweep ()
    { return neighbor_iterator (); };
    
    /// Return index of current octant across all
    /// trees on current process.
    p4est_locidx_t
    get_forest_oct_idx ()
    { return forest_oct_idx; };

    /// Return index of current octant in current
    /// tree on current process.
    p4est_locidx_t
    get_tree_oct_idx ()
    { return tree_oct_idx; };
    
    /// Return index of current octant across all
    /// trees on all processes.
    p4est_gloidx_t
    get_global_oct_idx ()
    {
      return forest_oct_idx +
        the_tmesh->p8est->global_first_quadrant[the_tmesh->rank];
    };
    
    /// Return index of current tree.
    p4est_locidx_t
    get_tree_idx ()
    { return tree_idx; };
    
    /// Update stored data.
    void
    update (p4est_topidx_t tree,
            p8est_quadrant_t *q);
    
    /// A pointer to the owning mesh is needed
    /// to get physical mapping.
    tmesh_3d                   *the_tmesh; 
    p8est_tree_t               *tree;
    p8est_quadrant_t           *the_octant;
    
    friend class tmesh_3d::octant_iterator;
    friend class tmesh_3d::neighbor_iterator;
    friend class tmesh_3d;
    
  private:

    sc_array_t           *toctants;
    p4est_locidx_t        num_octants;    // Q
    // Local and global indices for looping.
    p4est_topidx_t        tree_idx;         // tt
    p4est_locidx_t        forest_oct_idx;  // k 
    p4est_locidx_t        tree_oct_idx;    // q

    // True if current octant is a ghost.
    bool           is_ghost;
    // qtq index if current octant is a ghost, -1 otherwise.
    p4est_locidx_t qtq;
    
    /// Buffer used when quering coordinates. CESARE(da controllare)
    double                vxyz[24] = {0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0};//CESARE(24->P8EST_DIM*P8EST_CHILDREN)
    idx_t                 tbuff[8] = {0,0,0,0,0,0,0,0};//CESARE(8->P8EST_CHILDREN ?)
    bool                  hbuff[8] = {false,false,false,false,false,false,false,false};
    int                   pbuff[32] = {-1,-1, -1,-1, -1,-1, -1,-1,//CESARE(???)
                                       -1,-1, -1,-1, -1,-1, -1,-1,
                                       -1,-1, -1,-1, -1,-1, -1,-1,
                                       -1,-1, -1,-1, -1,-1, -1,-1};
  };

  /// Default constructor, set all pointers to nullptr.
  tmesh_3d (MPI_Comm _comm = MPI_COMM_WORLD)
    : p8est (nullptr), conn (nullptr),
      current_octant (this, 0, nullptr),
      lnodes (nullptr), mesh (nullptr), ghost(nullptr),
      mirror_data (nullptr), ghost_data (nullptr),
      comm (_comm), rank (0), size (1),
      replace_fun(userint_replace)
  {
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
  };

  /// Load a p8est and connectivity from a file.
  tmesh_3d (const char *filename, MPI_Comm _comm = MPI_COMM_WORLD)
    : tmesh_3d ()
  { load (filename); };

  /// Delete copy constructor.
  tmesh_3d (const tmesh_3d &) = delete;
  
  /// Delete assignment operator.
  tmesh_3d & operator=(const tmesh_3d &) = delete;
  
  /// Destructor.
  ~tmesh_3d ();

  /// Load a connectivity from a compressed binary
  /// Octave file then init the p8est.
  void
  read_connectivity (const char *filename,
                     int source = 0);

  /// Load a p8est and connectivity from a set of arrays
  /// then init the p8est.
  void
  read_connectivity (const double *p,
                     const p4est_topidx_t num_vertices,
                     const p4est_topidx_t *t,
                     const p4est_topidx_t num_trees,
                     int source = 0);

  /// Save the p8est and connectivity to a file.
  void
  save (const char *filename);

  /// Load the p8est and connectivity from a file.
  void
  load (const char *filename);

  /// Export exploded mesh to a vtk file for visualization.
  void
  vtk_export (const char *filename);

  /// Export nodal field f to a octbin.gz file for visualization.
  void
  octbin_export (const char * filename,
                 const std::vector<double> & f);
  
  /// Get an iterator to the first octant of the mesh.
  octant_iterator
  begin_octant_sweep ();

  /// Get a null octant iterator to signal end of the sweep.
  octant_iterator
  end_octant_sweep ()
  { return octant_iterator (); };

  /// Mark octant for refinement based on fun.
  void
  set_refine_marker
  (std::function<int (octant_iterator)> fun)
    {
      for (auto q = this->begin_octant_sweep ();
           q != this->end_octant_sweep ();
           ++q)
        if (fun (q))
          q->the_octant->p.user_int = std::abs(fun (q));
    };

  /// Mark octants for coarsening based on fun.
  void
  set_coarsen_marker
  (std::function<int (octant_iterator)> fun)
    {
      for (auto q = this->begin_octant_sweep ();
           q != this->end_octant_sweep ();
           ++q)
        if (fun (q))
          q->the_octant->p.user_int = -std::abs(fun (q));
    };
  
  /// Mark octants for refinement based on metrics.
  void set_metrics_marker (std::function<double (octant_iterator)>,
                           double, int max_depth = 5);

  /// Set functor to replace octants while being
  /// refined or coarsened.
  void
  set_replace_fun
  (std::function<std::vector<int> (std::vector<int>)> fun)
  { replace_fun = fun; };
  
  /// Refine marked octants, balance the octree and
  /// re-partition over the processors.
  void
  refine (int recursive = 0, int partforcoarsen = 1, int balance = 1);
  
  /// Refine marked octants based on metrics.
  void
  metrics_refine (idx_t max_elems = 0);

  /// Coarsen marked octants, balance the octree and
  /// re-partition over the processors.
  void
  coarsen (int recursive = 0, int partforcoarsen = 1, int balance = 1);

  /// Compute nodes numbering and octant neighbours.
  void
  update ();

  /// Send mirrors and receive ghosts from other ranks.
  void
  update_ghosts ();

  /// Return number of nodes owned by local process
  idx_t
  num_owned_nodes ()
  {
    if (! lnodes) update ();
    return lnodes->owned_count;
  };

  /// Return number of nodes of octants owned
  /// or shared by local process
  idx_t
  num_local_nodes ()
  {
      if (! lnodes) update ();
      return lnodes->num_local_nodes;
  };

  /// Return total number of octants owned by all process
  idx_t
  num_global_nodes ()
  {
    idx_t retval = 0;
    if (! lnodes) update ();
    for (int i = 0; i < size; ++i)
      retval += lnodes->global_owned_count[i];
    return retval;
  };

  /// Return number of octants owned by local process
  /// across all trees
  idx_t
  num_local_octants ()    
  {
      return p8est->local_num_quadrants;
  };

  /// Return number of octants owned by all processes
  /// across all trees
  idx_t
  num_global_octants ()    
  {
      return p8est->global_num_quadrants;
  };
  
  /// Replace fun based on octant user_int.
  static std::vector<int> userint_replace (std::vector<int>);
  
  /// P8EST pointers describing the tmesh,
  /// temporarily public until the API is stable.
  p8est_t              *p8est;
  p8est_connectivity_t *conn;
  octant_t              current_octant;
  p8est_lnodes_t       *lnodes;
  p8est_mesh_t         *mesh;
  p8est_ghost_t        *ghost;
  
  p4est_gloidx_t * mirror_data;
  p4est_gloidx_t * ghost_data;
  
  MPI_Comm comm;
  int      rank;
  int      size;
  
private:
  
  std::function<std::vector<int> (std::vector<int>)> replace_fun;

  static int
  refine_callback (p8est_t*, p4est_topidx_t, p8est_quadrant_t*);

  static int
  coarsen_callback (p8est_t*, p4est_topidx_t, p8est_quadrant_t* []);

  static void
  replace_callback (p8est_t*, p4est_topidx_t,
                    int, p8est_quadrant_t* [],
                    int, p8est_quadrant_t* []);
  
  int metrics_max_depth;
};

/// TODO : class for distributed arrays
template<class T, class data_location>
class
distributed_array
{

public :

  p4est_locidx_t
  num_local_entries ()
  { return data_location::num_local_entries (the_tmesh); }

  p4est_gloidx_t
  num_global_entries ()
  { return data_location::num_global_entries (the_tmesh); }

  T&
  index_local (p4est_locidx_t);

  const T&
  index_local (p4est_locidx_t) const;

  T&
  index_global (p4est_gloidx_t);

  const T&
  index_global (p4est_gloidx_t) const;

  void
  synchronize ();

  void
  all_reduce ();

private :

  tmesh_3d *the_tmesh;
  // should probably use some sort of smart pointer here
  T *data;
};


#endif /* TMESH_3D_H */

