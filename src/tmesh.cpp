/*
  Copyright (C) 2017 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file tmesh.cpp
  \brief Interface for p4est library
*/

#include <tmesh.h>
#include <array>


double
tmesh::quadrant_t::p (tmesh::idx_t ii, tmesh::idx_t jj) 
{
  double retval = vxyz[3*jj+ii];
  return (retval);
};

/** Decode the information from p{4,8}est_lnodes_t for a given element.
 *
 * \see p4est_lnodes.h for an in-depth discussion of the encoding.
 * \param [in] face_code         Bit code as defined in p{4,8}est_lnodes.h.
 * \param [out] hanging_corner   Undefined if no node is hanging.
 *                               If any node is hanging, this contains
 *                               one integer per corner, which is -1
 *                               for corners that are not hanging,
 *                               and the number of the non-hanging
 *                               corner on the hanging face/edge otherwise.
 *                               For faces in 3D, it is diagonally opposite.
 * \return true if any node is hanging, false otherwise.
 */
static const int    zero = 0;      /**< Constant zero. */
static const int    ones = P4EST_CHILDREN - 1;  /**< One bit per dimension. */
static int
lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN])
{
  if (face_code) {
    const int           c = (int) (face_code & ones);
    int                 i, h;
    int                 work = (int) (face_code >> P4EST_DIM);

    /* These two corners are never hanging by construction. */
    hanging_corner[c] = hanging_corner[c ^ ones] = -1;
    for (i = 0; i < P4EST_DIM; ++i) {
      /* Process face hanging corners. */
      h = c ^ (1 << i);
      hanging_corner[h ^ ones] = (work & 1) ? c : -1;
      work >>= 1;
    }
    return 1;
  }
  return 0;
}

void
tmesh::quadrant_t::update (p4est_topidx_t tree,
                           p4est_quadrant_t *q)
{
  p4est_quadrant_t node, parent;
  idx_t i;
  int hanging_corner[4];
  p4est_lnodes_t *ln = the_tmesh->lnodes;
  p4est_locidx_t lni;
  
  this->tree_idx = tree;
  this->the_quadrant = q;
  for (i = 0; i < 4; ++i)
    {
      p4est_quadrant_corner_node (this->the_quadrant, i, &node);
      p4est_qcoord_to_vertex (this->the_tmesh->conn, tree_idx,
                              node.x, node.y, &(vxyz[3 * i]));
    }

  if (ln != nullptr)
    {
      const p4est_locidx_t nloc = ln->num_local_nodes;
      
      for (i = 0; i < 4; ++i)
        {
          lni = ln->element_nodes[4 * forest_quad_idx + i];
          tbuff[i] = lni;
          hbuff[i] = false;
        }
            
      if (lnodes_decode2
          (ln->face_code[forest_quad_idx], hanging_corner))
        {
          for (i = 0; i < 4; ++i)
            if (hanging_corner[i] != -1)
              hbuff[i] = true;      
        }
    }
};

tmesh::idx_t
tmesh::quadrant_t::t (tmesh::idx_t i)
{ return tbuff[i]; };

bool
tmesh::quadrant_t::is_hanging (tmesh::idx_t i)
{ return hbuff[i]; };


tmesh::~tmesh ()
{
  p4est_destroy (this->p4est);
  p4est_connectivity_destroy (this->conn);
  if (!(this->lnodes == nullptr)) p4est_lnodes_destroy (this->lnodes);
};


/* Read a 2d p4est connectivity from a compressed octave
 * binary file. The file should contain a struct
 * named "msh" with fields "p" and "t". The former
 * should be the list of vertex coordinates while the
 * latter the list of element vertices.
 *
 * The following ordering is assumed for the nodes
 * of each quadrilateral element.
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
 */

static void
octbingz2connectivity (const char *filename, p4est_connectivity_t **conn)
{
    
  // load data from file
  octave_value tmp;
  octave_io_mode m = gz_read_mode;
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_load ("msh", tmp) == 0);

  Matrix p_matrix =
    tmp.scalar_map_value ().contents ("p").matrix_value ();

  Array<int> t_matrix =
    tmp.scalar_map_value ().contents ("t").array_value ();

  p4est_topidx_t num_vertices = p_matrix.cols (),
    num_trees = t_matrix.cols ();

  std::array<int, 4> v = {0, 0, 0, 0};
  int face = 0;

  *conn =
    p4est_connectivity_new (num_vertices, num_trees, 0, 0);

  for (int node = 1; node <= num_vertices; ++node)
    {
      (*conn)->vertices [3 * (node - 1) + 0] =
        p_matrix.fortran_vec ()[2 * (node - 1) + 0];
      
      (*conn)->vertices [3 * (node - 1) + 1] =
        p_matrix.fortran_vec ()[2 * (node - 1) + 1];
      
      (*conn)->vertices [3 * (node - 1) + 2] = 0;
    }

  for (int element_number = 1;
       element_number <= num_trees;
       ++element_number)
    {

      for (int n = 0; n < 4; ++n)
        v[n] = t_matrix.fortran_vec ()[5 * (element_number - 1) + n];

      (*conn)->tree_to_vertex[4 * (element_number - 1) + 0] = v[0] - 1;
      (*conn)->tree_to_vertex[4 * (element_number - 1) + 1] = v[1] - 1;
      (*conn)->tree_to_vertex[4 * (element_number - 1) + 2] = v[3] - 1;
      (*conn)->tree_to_vertex[4 * (element_number - 1) + 3] = v[2] - 1;

    }

  for (p4est_topidx_t tree = 0; tree < (*conn)->num_trees; ++tree)
    for (face = 0; face < 4; ++face)
      {
        (*conn)->tree_to_tree[4 * tree + face] = tree;
        (*conn)->tree_to_face[4 * tree + face] = face;
      }

  assert (p4est_connectivity_is_valid ((*conn)));
  p4est_connectivity_complete ((*conn));
};


void
tmesh::read_connectivity (const char *filename,
                          int source,
                          MPI_Comm comm)
{

  int rank, size;
  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &size);

  if (rank == source)
    octbingz2connectivity (filename, &conn);
  
  conn = p4est_connectivity_bcast (conn, source, comm);
  p4est = p4est_new (comm, conn, 0, NULL, NULL);
  p4est->user_pointer = this;

};

void
tmesh::save (const char *filename)
{ p4est_save (filename, p4est, 0); };

void
tmesh::load (const char *filename, MPI_Comm comm)
{ p4est = p4est_load (filename, comm, 0, 0, this, &conn); };

void
tmesh::vtk_export (const char *filename)
{
  p4est_vtk_context_t *context = p4est_vtk_context_new (p4est, filename);
  assert (context != nullptr);
  p4est_vtk_context_set_scale (context, 1.0);
  p4est_vtk_context_set_continuous (context, 1);
  context = p4est_vtk_write_header (context);
  context = p4est_vtk_write_cell_dataf (context, 1, 1, 1, 0, 0, 0, context);
  assert (p4est_vtk_write_footer (context) == 0);
};

void
tmesh::quadrant_iterator::reset ()
{
  p4est_t *p4 = data->the_tmesh->p4est;
  data->tree_idx         = p4->first_local_tree;
  data->forest_quad_idx  = 0;

  data->tree             = p4est_tree_array_index (p4->trees, data->tree_idx);
  data->tquadrants       = &(data->tree->quadrants);
  data->num_quadrants    = (p4est_locidx_t) data->tquadrants->elem_count;
  
  auto tmp = p4est_quadrant_array_index (data->tquadrants, data->forest_quad_idx);
  data->update (data->tree_idx, tmp);
  
};

tmesh::quadrant_iterator
tmesh::begin_quadrant_sweep ()
{
  quadrant_iterator qi (&current_quadrant);
  qi.reset ();
  return qi;
};

void 
tmesh::quadrant_iterator::operator++ ()
{

  p4est_t *p4 = data->the_tmesh->p4est;

  data->forest_quad_idx++;
  data->tree_quad_idx++;


  if (data->tree_quad_idx >= data->num_quadrants)
    {
      //std::cout << "tree_idx = " << data->the_tmesh->tree_idx << std::endl;
      data->tree_idx++;
      if ((data->tree_idx) > (p4->last_local_tree))
        {
          this->data = nullptr;
          return;
        }
      data->tree_quad_idx = 0;
      
      data->tree = p4est_tree_array_index (p4->trees, data->tree_idx);
      data->tquadrants = &(data->tree)->quadrants;

      data->num_quadrants = (p4est_locidx_t) data->tquadrants->elem_count;
    }

  auto tmp = p4est_quadrant_array_index (data->tquadrants, data->tree_quad_idx);
  data->update (data->tree_idx, tmp);
  
};

void
tmesh::refine (int recursive, int partforcoarsen)
{
  p4est_refine (p4est, recursive, refine_callback, nullptr);
  p4est_balance (p4est, P4EST_CONNECT_FULL, nullptr);
  p4est_partition (p4est, partforcoarsen, nullptr);

  if (! (lnodes == nullptr)) p4est_lnodes_destroy (lnodes);
  lnodes = nullptr;
}

void
tmesh::coarsen (int recursive, int partforcoarsen)
{
  p4est_coarsen (p4est, recursive, coarsen_callback, nullptr);
  p4est_balance (p4est, P4EST_CONNECT_FULL, nullptr);
  p4est_partition (p4est, partforcoarsen, nullptr);

  if (! (lnodes == nullptr)) p4est_lnodes_destroy (lnodes);
  lnodes = nullptr;  
};

void
tmesh::update ()
{
  auto ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  lnodes = p4est_lnodes_new (p4est, ghost, 1);
  p4est_ghost_destroy (ghost);
  ghost = nullptr;  
};

int
tmesh::refine_callback (p4est_t* p4, p4est_topidx_t tt, p4est_quadrant_t* qq)
{
  tmesh *tm = reinterpret_cast<tmesh*> (p4->user_pointer);
  tm->current_quadrant.update (tt, qq);
  quadrant_iterator qi (&(tm->current_quadrant));
  return tm->refine_marker (qi);
};

int
tmesh::coarsen_callback (p4est_t* p4, p4est_topidx_t tt, p4est_quadrant_t* qq[])
{
  tmesh *tm = reinterpret_cast<tmesh*> (p4->user_pointer);
  auto fun = [tm, tt, qq] (idx_t ii, quadrant_iterator& qi)
    { select_quad (tm, tt, qq, ii, qi); };
  return tm->coarsen_marker (fun);
};


void
tmesh::select_quad (tmesh *_tmesh,
                    p4est_topidx_t tree_idx,
                    p4est_quadrant_t* qt [],
                    idx_t ii, quadrant_iterator& qi)
{
  qi = quadrant_iterator (&(_tmesh->current_quadrant));
  qi->update (tree_idx, qt[ii]);
};
