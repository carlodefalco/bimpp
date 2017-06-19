/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


/* 
 * File:   tmesh.h
 * Author: Andrea Parini
 *
 * Created on 14 febbraio 2017, 9.25
 * 
 * 
 */

/*! \file tmesh.h
  \brief Interface for p4est library
*/

/*
 * tmesh is the p4est interface for the user: it contains all the necessary 
 * methods to allocate a p4est, refine or coarsen it, and access the mesh 
 * through quadrant iteration.
 * It also has methods for input/output in various formats.
 */

#ifndef TMESH_H
#define TMESH_H

#include <vector>
#include <array>
#include <mpi.h>
#include <functional>
#include <limits>
#include <tuple>

#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
#include <octave_file_io.h>
#include <cassert>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_vtk.h>
#include <p8est_mesh.h>
#include <p8est_extended.h>
#include <p8est_lnodes.h>
// 3D mode; must be included AFTER every other file
#include <p4est_to_p8est.h>

typedef p8est_connectivity_t* ConnectivityP;
typedef p8est_t* ForestP;
typedef p8est_tree_t * TreeP;
typedef p8est_quadrant_t * QuadrantP ;
typedef p8est_quadrant_t Quadrant ;
typedef p4est_topidx_t Int;
typedef p8est_lnodes_t* NodesP;


using namespace std;


//! Struct defining the state variables associated to each quadrant
typedef struct
{
  double eps;                 //!< Diffusion coefficient
  double psi[8];              //!< Scalar transport term
  double f[8];                //!< Forcing
  double u_sol[8];            //!< Solution
  Int tree;                   //!< Tree containing the quadrant
  vector< tuple<int,double> >bcs;  //!< Vector containing the boundary nodes of
                                  //!< the quadrant and their Dirichlet value
  vector< tuple<int,int,double> >dirichlet_bcs; //! Vector containing the 
                                  //!< assigned boundary conditions that are
                                  //!< this quadrant's competence (for faster
                                  //!< computation)
} user_data;


//! Struct used to collect quadrant and tree information
typedef struct
{
  QuadrantP quad;
  Int which_tree;
} tmesh_quadrant;

//! Struct used to define the present considered quadrant
/*!
 * Many of the tmesh class methods read and modify this structure, which holds
 * all the informations needed to manipulate the quadrants' informations. \n
 * Methods start_iterate and step_iterate are used for the quadrant iteration.
 */
typedef struct
{
  QuadrantP quad;      //!< pointer to the QuadrantP struct
  int index;           //!< index of the quadrant in the tree
  Int which_tree;      //!< index of the tree
  vector<int> gni;     //!< vector of 8 entries which keeps the global indexes 
                       //!< of the nodes
  int to_update_gni;   //!< 1 when a new quadrant is considered, 
                       //!< set to 0 by tmesh::t
  vector<int> lni;     //!< vector of 8 entries which keeps the processor local
                       //!< indexes of the nodes
  int to_update_lni;   //!< 1 when a new quadrant is considered, 
                       //!< set to 0 by tmesh::tloc
  vector< vector <int> > parents; //!< vector of 8 entries, 
                                  //!< for each entry it keeps 
                                  //!< the vector of the parents' global index 
                                  //!< for hanging nodes (or -1 if 
                                  //!< not hanging)
  vector< vector <int> > parents_loc; //!< vector of 8 entries, 
                                      //!< for each entry it keeps 
                                      //!< the vector of the parents' local
                                      //!< index for hanging nodes (or -1 if 
                                      //!< not hanging)
  
  vector< vector <vector<int> > > parents_borders; //!< For each node, keeps
                       //!< a vector of 4 elements for each possible parent node
                       //!< in which in stores the return of tmesh::e for that
                       //!< node
  
 int to_update_parents; //!< 1 when a new quadrant is considered, set to 0 by 
                         //!< tmesh::parent
  int to_update_parents_loc; //!< 1 when a new quadrant is considered, set to 0 by 
                             //!< tmesh::parent_loc
} current_quadrant;


//! tmesh class

class tmesh
{
private:
  
  int mpi_rank;          

  current_quadrant current;

  //! Callback function used to initialize the quadrants' user_data
  /*!
   * It needs to be passed to p8est_new to allocate the memory for the 
   * quadrants' user_data
   */
  static void
  p4_init
  (ForestP p8est, Int which_tree, QuadrantP quadrant);
  
  //! Callback function used for refinement
  /*!
   * It needs to be passed to p8est_refine, and evaluates marker_ref at each
   * quadrant of the forest
   */
  static int 
  callback_ref 
  (ForestP p8est, Int which_tree, QuadrantP q); 
  
  
  //! Callback function used for derefinement
  /*!
   * It needs to be passed to p8est_coarsen, and evaluates marker_deref at each
   * group of 8 sibling quadrants of the forest
   */
  static int 
  callback_deref
  (ForestP p8est, Int which_tree, QuadrantP q[]);
  
  
  //! Callback function used for refinement/derefinement
  /*!
   * It needs to be passed to p8est_refine_ext / p8est_coarsen_ext,
   * and allows to modify the parent quadrant's user_data before passing
   * it to its children, or vice versa. It calls the previously
   * set marker_replace
   */
  static void 
  callback_replace (ForestP p8est, 
                    Int which_tree,
                    int num_outgoing,
                    QuadrantP outgoing[],
                    int num_incoming,
                    QuadrantP incoming[]);
  
  
  
public:
  
  //! The MPI communicator
  MPI_Comm mpi_comm;
  
  //! The p4est connectivity field
  ConnectivityP conn;
  
  //! The p4est forest field
  ForestP forest;
  
  //! The p4est lnodes field
  NodesP nodes;
  
  //! Vector that associates to each tree its region number
  vector<int> tree_to_region;
  
  //! Region list
  vector<int> regions_idx;
  
  
  //! Marker for refinement
  /*!
   * Must be set by the user with set_refinement_marker, and it's evaluated
   * by callback_ref during refinement.
   * It must have the following structure:
   * \param [in] tmesh*     a valid pointer to the tmesh structure
   * \param [in] Int        the which_tree index
   * \param [in] QuadrantP  a pointer to the quadrant
   * \return 0 if the quadrant is not to refine or 1 if it is to refine
   */
  function <int(tmesh*, Int, QuadrantP)> marker_ref;
  
  
  //! Marker for derefinement
  /*!
   * Must be set by the user with set_derefinement_marker, and it's evaluated
   * by callback_deref during derefinement
   * It must have the following structure:
   * \param [in] tmesh*              a valid pointer to the tmesh structure
   * \param [in] Int                 the which_tree index
   * \param [in] array<QuadrantP,8>  the array of the sibling quadrants pointers
   * \return 0 if the quadrants are not to coarsen or 1 if they are to coarsen
   */
  function <int(tmesh*, Int, array <QuadrantP, 8> )> marker_deref;


  //! Marker for updating quadrants' user_data during refinement/derefinement
  /*!
   * Must be set by the user with set_replacement_marker, and it's evaluated
   * by callback_replace during refinement/derefinement
   * It must have the following structure:
   * \param [in] tmesh*              a valid pointer to the tmesh structure
   * \param [in] Int                 the which_tree index
   * \param [in] int                 the number of incoming quadrants
   * \param [in] vector<QuadrantP>   (num_incoming) the vector of incoming 
   *                                 quadrants pointers
   * \param [in] int                 the number of outgoing quadrants
   * \param [in] vector<QuadrantP>   (num_outgoing) the vector of outgoing 
   *                                 quadrants pointers
   */ 
  function <void(tmesh*, 
                 Int,
                 int, 
                 vector<QuadrantP>, 
                 int,
                 vector<QuadrantP>)> marker_replace;
    
  
  
  /* mesh creation/destruction methods */
  
  
  //! Basic constructor
  /*!
   * Initializes the tmesh class fields
   */
  tmesh
  (MPI_Comm mpi_comm_ = MPI_COMM_WORLD);
  
  
  //! Default destructor  
  ~tmesh() = default;
  
  
  //! Saves connectivity field and creates forest structure
  /*!
   * \param [in] conn_ the connectivity with which to create the forest
   */
  void
  build_forest(ConnectivityP conn_);
  
  
  //! Creates p4est lnodes structure
  /*!
   * This function must be called everytime after the forest is modified,
   * to update the lnodes structure
   */
  void
  build_p4nodes();
  
  
  //! Broadcast connectivity and tree_to_region and regions_idx fields
  /*! Broadcasts connectivity thanks to p4est_connectivity_bcast
   * and sends to other process the already allocated tree_to_region
   * and regions_idx vectors
   * \param [in] root  the processor that owns the connectivity 
   * \param [in] comm  a valid MPI communicator
   */
  void
  bcast_fields (int root, 
                MPI_Comm comm);
  
  
  //! Finalize mesh: should be called before MPI_Finalize
  void
  finalize();
  
  
  /* refinement/coarsening methods */
  
  
  
  //! Sets the refinement marker function
  /*!
   * \param [in] marker a valid lambda function with the prescribed signature
   */
  void 
  set_refinement_marker(function <int (tmesh*, Int, QuadrantP)> marker);
  
  
  //! Sets the derefinement marker function
  /*!
   * \param [in] marker a valid lambda function with the prescribed signature
   */
  void 
  set_derefinement_marker 
  (function <int (tmesh*, Int, array<QuadrantP, 8> )> marker);
  
  
  //! Sets the replacement marker function
  /*!
   * \param [in] marker a valid lambda function with the prescribed signature
   */
  void
  set_replace_marker 
  (function <void(tmesh*, 
                  Int, 
                  int, 
                  vector<QuadrantP>, 
                  int, 
                  vector<QuadrantP>)> marker);
  
  
  
  //! Refinement routine
  /*! 
   * Refines the forest based on the refinement marker, balances and
   * partitions the forest.
   * \param [in] to_balance         if 1 p8est_balance is called 
   *                                (keeps 2:1 ratio)
   * \param [in] recursive          if 1 refinement is recursive
   * \param [in] part_for_coarsen   if 1 quadrants families are never split
   *                                between processors during partition
   */
  void 
  refine(int to_balance,
         int recursive,
         int part_for_coarsen);
  
  
  //! Refinement routine with parameter interpolation
  /*! 
   * Refines the forest based on the refinement marker, interpolates user_data
   * parameters based on replace_marker, balances and partitions the forest.
   * \param [in] to_balance         if 1 p8est_balance is called
   * \param [in] recursive          if 1 refinement is recursive
   * \param [in] part_for_coarsen   if 1 quadrants families are never split
   *                                between processors during partition
   */
  void
  refine_interpolate(int to_balance, 
                     int recursive,
                     int part_for_coarsen);
  
  
  //! Derefinement routine
  /*! 
   * Derefines the forest based on the derefinement marker, balances and
   * partitions the forest.
   * \param [in] to_balance         if 1 p8est_balance is called
   *                                (keeps 2:1 ratio)
   * \param [in] recursive          if 1 derefinement is recursive
   * \param [in] part_for_coarsen   if 1 quadrants families are never split
   *                                between processors during partition
   */
  void 
  derefine(int to_balance,
           int recursive,
           int part_for_coarsen);
  
  
  //! Derefinement routine with parameter interpolation
  /*! 
   * Derefines the forest based on the derefinement marker interpolates 
   * user_data parameters based on replace_marker, 
   * balances and partitions the forest.
   * \param [in] to_balance         if 1 p8est_balance is called
   *                                (keeps 2:1 ratio)
   * \param [in] recursive          if 1 derefinement is recursive
   * \param [in] part_for_coarsen   if 1 quadrants families are never split
   *                                between processors during partition
   */
  void 
  derefine_interpolate(int to_balance,
                       int recursive,
                       int part_for_coarsen);
  
  
  /* mesh save/load methods*/
  
  
  //! Reads a connectivity from a compressed octave binary file
  /*!
   * Reads p and t mesh fields from the passed file name and creates a
   * connectivity
   * \param [in] filename the c-style string name of the file
   * \return a pointer to a valid connectivity
   */
  ConnectivityP 
  read_connectivity_oct_bin_gz(const char* filename);
  
  
  
  //! Reads a connectivity from a compressed octave binary file
  /*!
   * Reads p and t mesh fields from the passed filename input and creates a
   * connectivity
   * \param [in] filename the string name of the file
   * \return a pointer to a valid connectivity
   */
  ConnectivityP 
  read_connectivity_oct_bin_gz(string filename);
  
  
  //! Writes the passed p and t vectors to a compressed octave binary file
  /*!
   * Takes the vertices vector p and the element vector t and writes them 
   * to a compressed octave binary file named filename.
   * The following ordering is assumed for the nodes
   * of each element. This is the ordering for the nodes for the lesser z
   * coordinate. The nodes with the greater z coordinates are ordered in the 
   * same manner after these ones.
   *
   *   4                   3
   *   +-------------------+
   *   |                   |
   *   |                   |
   *   |                   |
   *   |                   |
   *   |                   |
   *   |                   |
   *   +-------------------+
   *   1                   2
   * \param [in] filename the c-style string name of the file
   * \param [in] p        (3*num_vertices) the point vector
   * \param [in] t        (8*num_elems) the element vector
   * \return              0 if writing was successful
   */
  bool 
  write_connectivity_oct_bin_gz(const char* filename, 
                                vector<double> p, 
                                vector<int> t);
  
  
  //! Writes the passed p and t vectors to a compressed octave binary file
  /*!
   * Takes the vertices vector p and the element vector t and writes them 
   * to a compressed octave binary file named filename.
   * The following ordering is assumed for the nodes
   * of each element. This is the ordering for the nodes for the lesser z
   * coordinate. The nodes with the greater z coordinates are ordered in the 
   * same manner after these ones.
   *
   *   4                   3
   *   +-------------------+
   *   |                   |
   *   |                   |
   *   |                   |
   *   |                   |
   *   |                   |
   *   |                   |
   *   +-------------------+
   *   1                   2
   * \param [in] filename the string name of the file
   * \param [in] p        (3*num_vertices) the point vector
   * \param [in] t        (8*num_elems) the element vector. See above for node 
   *                      numbering convention that must be given in the 
   *                      t vector
   * \return              0 if writing was successful
   */
  bool 
  write_connectivity_oct_bin_gz(string filename,
                                vector<double> p, 
                                vector<int> t);
  
  
  //! Creates a vtk file with the forest information
  /*!
   * \param [in] file_out  the string name of the file
   */
  void
  save_vtk(string file_out);
  
  
  
  
  //! Saves the mesh and the function f on a compressed binary octave file
  /*!
   * \param[in] f_vec     a vector (8*num_local_quadrants) with the nodal values
   *                               of a function to be saved
   * \return 0  if the saving was successful
   */
  int
  save_oct_bin (string filename, 
                vector<double> f_vec); 
  
  
  /* mesh access methods */
  
  //! Returns the number of owned quadrants
  /*!
   * \return  the processor local owned number of quadrants
   */
  int
  num_owned_quadrants();
  
  
  //! Returns the number of owned nodes
  /*!
   * \return  the processor local number of owned nodes (does not count shared
   *          nor hanging nodes)
   */
  int
  num_owned_nodes();
  
  
  //! Returns the number of local nodes, counting shared ones
  /*!
   * \return  the processor local number of local nodes (does count shared ones
   *          but not hanging nodes)
   */
  int 
  num_local_nodes();
  
  
  //! Returns the total number quadrants on all processors
  /*!
   * \return  the number of quadrants on all processors
   */
  int
  nelem();
  
  
  //! Returns the total number of independent nodes on all processors
  /*!
   * \return  the number of independent nodes on all processors
   */
  int
  nnodes();
  
  
  //! Starts the quadrant iteration of the tree indexed by which_tree
  /*!
   * Starts the quadrant iteration on the which_tree tree by setting the
   * current quadrant information. current.quad is set to the first quadrant 
   * of the tree
   * \param [in] which_tree  the tree on which the iteration is started 
   */
  void
  start_iterate(Int which_tree);
  
  
  //! Advances the quadrant iteration to the next quadrant
  /*!
   * Advances the quadrant iteration of the which_tree tree to the next 
   * quadrant, must be called after start_iterate.
   * If called after the end of a tree, it will just return
   */
  void
  step_iterate();
  
  
  //! Get method for the current quadrant
  /*!
   * Saves the current_quadrant information to a current_quadrant struct
   * \return  a valid current_quadrant struct
   */
  current_quadrant
  get_current_information();
  
  
  //! Returns the global index of the inode-th vertex of the current quadrant
  /*!
   * Saves in current.gni the global node indexes of the
   * 8 vertices of the current quadrant if it has never been called on the
   * current quadrant, then it returns the inode-th element of the saved
   * current.gni. MUST be called after start_iterate, and step_iterate must
   * be used to advance to the next quadrant.
   * (Note that if called on hanging nodes it returns a non sense value, as
   * hanging nodes do not have a global node index)
   * \param [in] inode  the number 0-7 of the wanted node index
   * \return  the global node index of the vertex
   */
  int
  t (const int inode);  
  
  
  //! Returns the local index of the inode-th vertex of the quadrant
  //! with local index qidx
  /*!
   * Saves in current.lni the local node indexes of the
   * 8 vertices of the current quadrant if it has never been called on the
   * current quadrant, then it returns the inode-th element of the saved
   * current.lni. MUST be called after start_iterate, and step_iterate must
   * be used to advance to the next quadrant.
   * \param [in] inode  the number 0-7 of the wanted node index.
   * \return  the local node index of the vertex
   */
  int
  tloc (const int inode);

  
  //! Returns the coordinates of the inode-th vertex of the current quadrant
  /*!
   * \param [in] idir   the wanted coordinate 0-2
   * \param [in] inode  the number 0-7 of the wanted node coordinate
   * \return  the coordinate of the node 
   */
  double
  p (const int idir, const int inode);

  
  //! Returns the coordinates of the vertex with global index inode
  /*!
   * \param [in] direction   the wanted coordinate 0-2
   * \param [in] inode       the global node index of the vertex
   * \return  the coordinate of the node
   */
  /* used for output only*/
  double
  p_out (const int direction, const int inode);

  
  //! Returns the global index of the inode-th locally owned vertex
  /*!
   * Must not be called on hanging nodes
   * \param [in] inode  the local node index of the node
   * \return  the global node index of the node
   */
  int
  l2g (const int inode);

  
  //! Returns the list of indexes (0-5) of the connectivity faces on which 
  //! the inode-th vertex of a quadrant is lying,
  //! return -1 if the vertex is not on an edge of the connectivity
  /*!
   * Returns a vector whose entries are numbered between 0 and 5, meaning that
   * the vertex inode of the tmesh_quadrant tq lies on the connectivity
   * faces of the tree containing the quadrant. Returns -1 if the
   * node lies on none of the tree faces. The method does not require an
   * active iteration, and can be called on any valid tmesh_quadrant.
   * \param [in] inode  the inode-th vertex of the quadrant
   * \param [in] tq     a struct containing the quadrant and its tree
   * \return  a vector containing entries between 0-5 or -1
   */
  vector<int>
  e (const int inode, tmesh_quadrant tq);
  
  
  //! Returns the index of the region to which the tree belongs
  /*!
   * \param [in] tq   a struct containing the quadrant and its tree
   * \return  the index of the region containing the tree
   */
  int
  r (tmesh_quadrant tq);
  
  

  //! Return the global index of the current quadrant on its processor
  /*!
   * Returns the index of the current quadrant between all 
   * quadrants on the processor.
   * Must be called on an active iteration, using start_iterate and 
   * step_iterate.
   * \return  the index of the current quadrant across all processor owned 
   *          quadrants
   */
  int
  idx ();
  
  
  //! Return the global index of the current quadrant across all processors
  /*!
   * Returns the index of the current quadrant between all 
   * quadrants on all processors.
   * Must be called on an active iteration, using start_iterate and 
   * step_iterate.
   * \return  the index of the current quadrant across all quadrants
   */
  int
  glo_idx ();
  
  
  
  

  //! Returns the global node index of a parent of a hanging node
  /*!
   * Saves in current.parents the global node indexes of the
   * parents of the 8 vertices of the current quadrant, 
   * if it has never been called on the current quadrant, then it returns 
   * the i_parent-th parent of the quad_j_node-th vertex of the current 
   * quadrant.
   * MUST be called after start_iterate, and step_iterate must
   * be used to advance to the next quadrant.
   * \param [in] i_parent            the number 0-4 of the wanted node parent
   * \param [in] i_node              the index 0-7 of the node
   * \return  the global node index of the parent if hanging,
   *          or -1 if the node is not hanging,
   *          or -2 if i_parent >= 4
   */
  int
  parent (const int i_parent,
          const int i_node);
  
  
  
  //! Returns the local node index of a parent of a hanging node
  /*!
   * Saves in current.parents the local node indexes of the
   * parents of the 8 vertices of the current quadrant, 
   * if it has never been called on the current quadrant, then it returns 
   * the i_parent-th parent of the quad_j_node-th vertex of the current 
   * quadrant.
   * MUST be called after start_iterate, and step_iterate must
   * be used to advance to the next quadrant.
   * \param [in] i_parent            the number 0-4 of the wanted node parent
   * \param [in] i_node              the index 0-7 of the node
   * \return  the local node index of the parent if hanging,
   *          or -1 if the node is not hanging,
   *          or -2 if i_parent >= 4
   */
  int
  parent_loc (const int i_parent,
              const int i_node);

  
  //! Decode the face_code into hanging face information.
  /*
   * \see p4est_lnodes.h for an in-depth discussion of the encoding.
   * \param[in] face_code as in the p8est_lnodes_t structure.
   * \param[out] hanging_face: if there are hanging faces or edges,
   *             hanging_face = -1 if the face is not hanging,
   *                          = the corner of the full face that it touches:
   *                            e.g. if face = i and hanging_face[i] =
   *                            j, then the interpolation operator corresponding
   *                            to corner j should be used for that face.
   *             note: not touched if there are no hanging faces or edges.
   * \param[out] hanging_edge: if there are hanging faces or edges,
   *             hanging_edge = -1 if the edge is not hanging,
   *                          =  0 if the edge is the first half of a full edge,
   *                              but neither of the two faces touching the
   *                              edge is hanging,
   *                          = 1 if the edge is the second half of a full edge,
   *                              but neither of the two faces touching the
   *                              edge is hanging,
   *                          = 2 if the edge is the first half of a full edge
   *                              and is on the boundary of a full face,
   *                          = 3 if the edge is the second half of a full edge
   *                              and is on the boundary of a full face,
   *                          = 4 if the edge is in the middle of a full face.
   *             note: not touched if there are no hanging faces or edges.
   * \return             true if any face or edge is hanging, false otherwise.
   */
  int
  lnodes_decode 
  (p8est_lnodes_code_t face_code, int hanging_face[6], int hanging_edge[12]);
 
  
  //! Decode the face_code into hanging face information.
  /*
   * \see p4est_lnodes.h for an in-depth discussion of the encoding
   * \param [in] face_code as in the p8est_lnodes_t structure.
   * \param [out] hanging_corner   Undefined if no node is hanging.
   *                               If any node is hanging, this contains
   *                               one integer per corner, which is -1
   *                               for corners that are not hanging,
   *                               and the number of the non-hanging
   *                               corner on the hanging face/edge otherwise.
   *                               For faces in 3D, it is diagonally opposite.
   * \return             true if any face or edge is hanging, false otherwise.
   */
  int
  lnodes_decode2
  (p8est_lnodes_code_t face_code, int hanging_corner[P4EST_CHILDREN]);
  
  
  /* utility methods */
  
  //! Method for safe double confrontation
  /*!
   * Returns true if doubles a and b are almost equal, where almost means that
   * the absolute value of their difference is smaller than a number that is
   * proportional to the order of magnitude of their sum. The method is needed
   * for the cases in which confrontation between values processed by the
   * library is required, and these values might be slightly different from 
   * the expected ones, so confronting them using equal operator is unsafe and
   * leads to errors
   * \param [in] a,b  doubles to be confronted
   * \return    true if they are equal 
   */
  int
  almost_equal(const double a, const double b);
  
  
  
};


#endif /* TMESH_H */

