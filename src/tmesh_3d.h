
/*! \file tmesh_3d.h
  \brief Interface for p4est library
*/

#ifndef TMESH_3D_H
#define TMESH_3D_H

#include <mpi.h>

#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_vtk.h>

#include <bim_distributed_vector.h>
#include <bim_ordering.h>

#include <functional>
#include <array>
#include <utility>
#include <vector>

/// C++ interface class for p4est 3d quadrant meshes.
class
tmesh_3d
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
    { return ! ((*this) == other); };

    /// Default constructor.
    quadrant_iterator (quadrant_t *_data = nullptr) :
      data (_data)
    { };

    /// Move to first forest quadrant
    void
    reset ();

  protected:
    quadrant_t *data;
  };

  /// Iterator to sweep through the quadrants of a tmesh.
  /// This is essentially a decorator of quadrant_t*.
  class
  neighbor_iterator : public quadrant_iterator
  {
  public:
    /// Get next neighbor.
    void
    operator++ ();

    /// Default constructor.
    neighbor_iterator (quadrant_t *_data = nullptr,
		       int _face_idx = -1) :
      quadrant_iterator (_data),
      face_neighbor (new p8est_mesh_face_neighbor_t),
      face_idx (_face_idx)
    { };

    /// Destructor.
    ~neighbor_iterator ()
    {
      if (data != nullptr)
	delete data;

      if (face_neighbor != nullptr)
	delete face_neighbor;
    }

    /// Get the face index associated to the current neighbor.
    int
    get_face_idx ()
    { return face_idx; };

    friend class tmesh_3d::quadrant_t;

  private:
    p8est_mesh_face_neighbor_t * face_neighbor; // mfn

    int face_idx; /// Face index in 0...5 (-1 if not defined).
  };

  /// C++ interface class to access properties of the
  /// current quadrant.
  class
  quadrant_t
  {

  public:

    /// Simple constructor needs at least a pointer
    /// to the container tmesh.
    quadrant_t (tmesh_3d *_tmesh,
		p4est_topidx_t _tree = 0,
		p8est_quadrant_t *_quadrant = nullptr) :
      the_tmesh (_tmesh), the_quadrant (_quadrant), tree_idx (_tree),
      is_ghost (false), qtq (-1)
    { };

    /// Get the i-th coordinate of the j-th vertex.
    double
    p (idx_t i, idx_t j);

    /// Get the i-th coordinate of the centroid.
    double
    centroid (idx_t i);

    /// Get the i-th coordinate of the centroid of the j-th face
    double
    face_centroid (idx_t i, int j);

    /// Get rank-local index of the i-th vertex, or global index for ghosts.
    idx_t
    t (idx_t i);

    /// Get global index of the i-th vertex
    idx_t
    gt (idx_t i);

    /// True if the i-th vertex is hanging.
    bool
    is_hanging (idx_t i);

    /// Get number of parents of the i-th vertex (0 if not hanging).
    int
    num_parents (idx_t i);

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
    /// Index of the face of the current tree
    //  on which the i-th vertex lies, return
    //  NOT_ON_BOUNDARY if an interior vertex.
    idx_t
    e (idx_t i);

    /// face_nodes[ii] is the list of vertices
   //  on the ii-th face of a quadrant,in local
   //  numbering.
   const std::array<std::vector<idx_t>, 7> 
   face_nodes  {std::vector<idx_t>{0, 2, 4, 6}, 
                 std::vector<idx_t>{1, 3, 5, 7}, 
                 std::vector<idx_t>{0, 1, 4, 5}, 
                 std::vector<idx_t>{2, 3, 6, 7}, 
                 std::vector<idx_t>{0, 1, 2, 3}, 
                 std::vector<idx_t>{4, 5, 6, 7},
                 std::vector<idx_t>{}};

    /// List of the nodes of the current quadrant
    //  that lie on the i-th face of the tree.
    const std::vector<idx_t> &
    ef (idx_t i);

    /// Get an iterator to the first neighbor
    /// of the current quadrant.
    neighbor_iterator
    begin_neighbor_sweep ();

    /// Get a null quadrant iterator to signal end of the sweep.
    neighbor_iterator
    end_neighbor_sweep ()
    { return neighbor_iterator (); };

    /// Return index of current quadrant across all
    /// trees on current process.
    p4est_locidx_t
    get_forest_quad_idx ()
    { return forest_quad_idx; };

    /// Return index of current quadrant in current
    /// tree on current process.
    p4est_locidx_t
    get_tree_quad_idx ()
    { return tree_quad_idx; };

    /// Return index of current quadrant across all
    /// trees on all processes.
    p4est_gloidx_t
    get_global_quad_idx ()
    {
      return forest_quad_idx +
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
    p8est_quadrant_t           *the_quadrant;

    friend class tmesh_3d::quadrant_iterator;
    friend class tmesh_3d::neighbor_iterator;
    friend class tmesh_3d;

  private:

    sc_array_t           *tquadrants;
    p4est_locidx_t        num_quadrants;    // Q
    // Local and global indices for looping.
    p4est_topidx_t        tree_idx;         // tt
    p4est_locidx_t        forest_quad_idx;  // k
    p4est_locidx_t        tree_quad_idx;    // q

    // True if current quadrant is a ghost.
    bool           is_ghost;
    // qtq index if current quadrant is a ghost, -1 otherwise.
    p4est_locidx_t qtq;

    /// Buffer used when quering coordinates.
    double vxyz [3 * 8] = {0,0,0, 0,0,0, 0,0,0, 0,0,0,
			   0,0,0, 0,0,0, 0,0,0, 0,0,0};

    /// Buffer for index of the i-th vertex.
    idx_t  tbuff[8]     = {0,0,0,0,0,0,0,0};
    /// Buffer for num. of parents of hanging nodes(0 if not hanging).
    int    hbuff[8]     = {0,0,0,0,0,0,0,0};
    /// Buffer for parents' t(-1 if not hanging or less than 4 parents).
    int    pbuff[4 * 8] = {-1,-1,-1,-1, -1,-1,-1,-1,
			   -1,-1,-1,-1, -1,-1,-1,-1,
			   -1,-1,-1,-1, -1,-1,-1,-1,
			   -1,-1,-1,-1, -1,-1,-1,-1};
  };

  /// Struct for p8est user_data.
  struct data_t
  {
    /// Number of refinement steps to be performed.
    /// A negative number is used to mark for coarsening.
    int refine_count;

    /// Interpolation indices, i.e. the indices
    /// associated to interp_coeff columns.
    // std::array<tmesh_3d::idx_t, 8> interp_idx;

    /// Interpolation coefficients at the eight vertices.
    //  std::array<std::array<double, 8>, 8> interp_coeff;
  };

  /// Default constructor, set all pointers to nullptr.
  tmesh_3d (MPI_Comm _comm = MPI_COMM_WORLD)
    : p8est (nullptr), conn (nullptr),
      current_quadrant (this, 0, nullptr),
      lnodes (nullptr), mesh (nullptr), ghost (nullptr),
      mirror_data (nullptr), ghost_data (nullptr),
      comm (_comm), rank (0), size (1),
      replace_fun (user_data_replace)
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
  tmesh_3d &
  operator= (const tmesh_3d &) = delete;

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
		 const std::vector<double> & f,
		 const ordering& ord = default_ord);

  /// Export nodal field f to a octbin.gz file for visualization.
  void
  octbin_export (const char * filename,
		 const distributed_vector & f,
		 const ordering& ord = default_ord);

  /// Export quadrant field f to a octbin.gz file for visualization.
  void
  octbin_export_quadrant (const char * filename,
			  const std::vector<double> & f);

  /// Get an iterator to the first quadrant of the mesh.
  quadrant_iterator
  begin_quadrant_sweep ();

  /// Get a null quadrant iterator to signal end of the sweep.
  quadrant_iterator
  end_quadrant_sweep ()
  { return quadrant_iterator (); };

  /// Mark quadrants for refinement based on fun.
  void
  set_refine_marker
  (std::function<int (quadrant_iterator)> fun)
  {
    tmesh_3d::data_t * data;
    int val = 0;

    for (auto q = this->begin_quadrant_sweep ();
	 q != this->end_quadrant_sweep ();
	 ++q)
      {
	// set_interpolation_matrix (q);

	val = fun (q);
	if (val)
	  {
	    data = static_cast<tmesh_3d::data_t *> (q->the_quadrant->p.user_data);
	    data->refine_count = std::abs (val);
	  }
      }
  };

  /// Mark quadrants for coarsening based on fun.
  void
  set_coarsen_marker
  (std::function<int (quadrant_iterator)> fun)
  {
    tmesh_3d::data_t * data;
    int val = 0;

    for (auto q = this->begin_quadrant_sweep ();
	 q != this->end_quadrant_sweep ();
	 ++q)
      {
	// set_interpolation_matrix (q);

	val = fun (q);
	if (val)
	  {
	    data = static_cast<tmesh_3d::data_t *> (q->the_quadrant->p.user_data);
	    data->refine_count = -std::abs (val);
	  }
      }
  };

  /// Mark quadrants for refinement based on metrics.
  void
  set_metrics_marker (std::function<double (quadrant_iterator)>,
		      double, int max_depth = 5,
		      int n_refine = 0,
		      int n_coarsen = 0);

  /// Set functor to replace quadrants while being
  /// refined or coarsened.
  void
  set_replace_fun
  (std::function<void (std::vector<tmesh_3d::data_t *>, std::vector<tmesh_3d::data_t>&)> fun)
  { replace_fun = fun; };

  /// Refine marked quadrants, balance the octree and
  /// re-partition over the processors.
  void
  refine (int recursive = 0, int partforcoarsen = 1, int balance = 1);

  /// Refine marked quadrants based on metrics.
  void
  metrics_refine (idx_t max_elems = 0);

  /// Coarsen marked quadrants, balance the octree and
  /// re-partition over the processors.
  void
  coarsen (int recursive = 0, int partforcoarsen = 1, int balance = 1);

  /// Compute nodes numbering and quadrant neighbours.
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

  /// Return number of nodes of quadrants owned
  /// or shared by local process
  idx_t
  num_local_nodes ()
  {
    if (! lnodes) update ();
    return lnodes->num_local_nodes;
  };

  /// Return total number of quadrants owned by all process
  idx_t
  num_global_nodes ()
  {
    idx_t retval = 0;
    if (! lnodes) update ();
    for (int i = 0; i < size; ++i)
      retval += lnodes->global_owned_count[i];
    return retval;
  };

  /// Return number of quadrants owned by local process
  /// across all trees
  idx_t
  num_local_quadrants ()
  {
    return p8est->local_num_quadrants;
  };

  /// Return number of quadrants owned by all processes
  /// across all trees
  idx_t
  num_global_quadrants ()
  {
    return p8est->global_num_quadrants;
  };

  /// Replace fun based on quadrant user_int.
  static std::vector<int>
  userint_replace (std::vector<int>);

  /// Replace fun based on quadrant user_data.
  static void
  user_data_replace (std::vector<tmesh_3d::data_t *>, std::vector<tmesh_3d::data_t>&);

  /// P8EST pointers describing the tmesh,
  /// temporarily public until the API is stable.
  p8est_t              *p8est;
  p8est_connectivity_t *conn;
  quadrant_t            current_quadrant;
  p8est_lnodes_t       *lnodes;
  p8est_mesh_t         *mesh;
  p8est_ghost_t        *ghost;

  p4est_gloidx_t       *mirror_data;
  p4est_gloidx_t       *ghost_data;

  MPI_Comm comm;
  int      rank;
  int      size;

private:
  std::function<void (std::vector<tmesh_3d::data_t *>,
		      std::vector<tmesh_3d::data_t>&)> replace_fun;

  static int
  refine_callback (p8est_t*, p4est_topidx_t, p8est_quadrant_t*);

  static int
  coarsen_callback (p8est_t*, p4est_topidx_t, p8est_quadrant_t* []);

  static void
  replace_callback (p8est_t*, p4est_topidx_t,
		    int, p8est_quadrant_t* [],
		    int, p8est_quadrant_t* []);

  // void
  // set_interpolation_matrix (tmesh_3d::quadrant_iterator &);

  int metrics_max_depth;
};


// void
// make_connectivity_3d (const p4est_topidx_t num_trees[3], const double step[3],
// 		      double *& p, p4est_topidx_t & num_vertices,
// 		      p4est_topidx_t *& t, p4est_topidx_t & total_num_trees,
// 		      std::vector<std::pair<p4est_topidx_t, p4est_topidx_t>> & bcells);
void
make_connectivity_3d (const p4est_topidx_t num_trees[3],
          const double step[3],
          std::unique_ptr<double[]> & p,
          p4est_topidx_t & num_vertices,
          std::unique_ptr<p4est_topidx_t[]> & t,
          p4est_topidx_t & total_num_trees,
          std::vector<std::pair<p4est_topidx_t, p4est_topidx_t>> &bcells);

#endif /* TMESH_3D_H */

