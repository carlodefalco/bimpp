/*! \file tmesh_3d.cpp
  \brief Interface for p8est library
*/

#include <tmesh_3d.h>
#include <array>

double
tmesh_3d::octant_t::p (tmesh_3d::idx_t ii, tmesh_3d::idx_t jj) 
{
  double retval = vxyz[3*jj+ii];
  return (retval);
};

double
tmesh_3d::octant_t::centroid (tmesh_3d::idx_t ii) 
{
  double retval = 0;
  
  for (tmesh_3d::idx_t c = 0; c < P8EST_CHILDREN; ++c)
    retval += this->p(ii, c);
  retval /= P8EST_CHILDREN;
  
  return (retval);
};

double
tmesh_3d::octant_t::face_centroid (tmesh_3d::idx_t ii, int jj)
{
  double retval = 0;
  
  for (tmesh_3d::idx_t c = 0; c < P8EST_HALF; ++c)
    retval += this->p(ii, p8est_face_corners[jj][c]);
  retval /= P8EST_HALF;
  
  return (retval);
}

/** Decode the information from p8est_lnodes_t for a given element.
 *
 * \see p8est_lnodes.h for an in-depth discussion of the encoding.
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
static const int  zero = 0;      /**< Constant zero. */
static const int  ones = P8EST_CHILDREN - 1;  /**< One bit per dimension. */
static const int *corner_to_hanging[P8EST_CHILDREN];
static const int  corner_num_hanging[P8EST_CHILDREN] = { 1, 2, 2, 4,
                                                         2, 4, 4, 1 };

static int
lnodes_decode2 (p8est_lnodes_code_t face_code,
                int hanging_corner[P8EST_CHILDREN])
{
  if (face_code) {
    const int           c = (int) (face_code & ones);
    int                 i, h;
    int                 work = (int) (face_code >> P8EST_DIM);

    /* These two corners are never hanging by construction. */
    hanging_corner[c] = hanging_corner[c ^ ones] = -1;
    for (i = 0; i < P8EST_DIM; ++i) {
      /* Process face hanging corners. */
      h = c ^ (1 << i);
      hanging_corner[h ^ ones] = (work & 1) ? c : -1;
      /* Process edge hanging corners. */
      hanging_corner[h] = (work & P8EST_CHILDREN) ? c : -1;
      work >>= 1;
    }
    return 1;
  }
  return 0;
}

void 
tmesh_3d::octant_iterator::operator++ ()
{

  p8est_t *p8 = data->the_tmesh->p8est;

  data->forest_oct_idx++;
  data->tree_oct_idx++;

  if (data->tree_oct_idx >= data->num_octants)
    {
      data->tree_idx++;
      
      if ((data->tree_idx) > (p8->last_local_tree))
        {
          this->data = nullptr;
          return;
        }
      data->tree_oct_idx = 0;
      
      data->tree = p8est_tree_array_index (p8->trees, data->tree_idx);
      data->toctants = &(data->tree)->quadrants;

      data->num_octants =
        (p4est_locidx_t) data->toctants->elem_count;
    }

  auto tmp = p8est_quadrant_array_index (data->toctants,
                                         data->tree_oct_idx);
  data->update (data->tree_idx, tmp);
  
};

void 
tmesh_3d::neighbor_iterator::operator++ ()
{
  p4est_topidx_t which_tree;
  p4est_locidx_t which_oct;
  int nface, nrank;
  tmesh_3d *tmsh = data->the_tmesh;
  p8est_t *p8 = tmsh->p8est;
  
  p8est_quadrant_t * neighbor =
    p8est_mesh_face_neighbor_next (face_neighbor, &which_tree,
                                   &which_oct, &nface, &nrank);
  
  if (neighbor != nullptr)
    {
      p8est_tree_t * tree =
        p8est_tree_array_index (p8->trees,
                                which_tree);
      
      // If non-ghost.
      if (face_neighbor->current_qtq <
          tmsh->num_local_octants ())
        {
          data->is_ghost = false;
          data->qtq = -1;
          
          data->forest_oct_idx = tree->quadrants_offset + which_oct;
          data->tree_oct_idx = which_oct;
        }
      // If ghost.
      else
        {
          data->is_ghost = true;
          data->qtq = face_neighbor->current_qtq;
          
          data->forest_oct_idx =
            neighbor->p.piggy3.local_num +
            (p8->global_first_quadrant[nrank] -
             p8->global_first_quadrant[tmsh->rank]);
          
          data->tree_oct_idx = data->forest_oct_idx -
            tree->quadrants_offset;
        }
      
      this->face_idx = nface;
      
      data->update (which_tree, neighbor);
    }
  else
    {
      data = nullptr;
      this->face_idx = -1;
    }
};

void
tmesh_3d::octant_t::update (p4est_topidx_t tree,
                           p8est_quadrant_t *q)
{
  p8est_quadrant_t node, parent;
  idx_t i, j;
  int hanging_corner[8];
  p8est_lnodes_t *ln = the_tmesh->lnodes;
  p4est_locidx_t lni;

  int c, h, num_parents;
  const int *base_corner;
  
  corner_to_hanging[0]        = &zero;
  corner_to_hanging[1]        = p8est_edge_corners[0];
  corner_to_hanging[2]        = p8est_edge_corners[4];
  corner_to_hanging[3]        = p8est_face_corners[4];
  corner_to_hanging[4]        = p8est_edge_corners[8];
  corner_to_hanging[ones - 2] = p8est_face_corners[2];
  corner_to_hanging[ones - 1] = p8est_face_corners[0];
  corner_to_hanging[ones]     = &ones;
  
  this->tree_idx = tree;
  this->the_octant = q;
  
  for (i = 0; i < 8; ++i)
    {
      p8est_quadrant_corner_node (this->the_octant, i, &node);
      p8est_qcoord_to_vertex (this->the_tmesh->conn, tree_idx,
                              node.x, node.y, node.z, &(vxyz[3 * i]));
    }

  if (ln != nullptr)
    {
      // Non-ghost elements.
      if (! is_ghost)
        {
          for (i = 0; i < 8; ++i)
            {
              tbuff[i] = ln->element_nodes[8 * forest_oct_idx + i];
              hbuff[i] = false;
              pbuff[4 * i] = -1;
              pbuff[4 * i + 1] = -1;
              pbuff[4 * i + 2] = -1;
              pbuff[4 * i + 3] = -1;
            }

          bool any_hanging =
            lnodes_decode2 (ln->face_code[forest_oct_idx],
                            hanging_corner);
          if (any_hanging)
            for (i = 0; i < 8; ++i)
              if (hanging_corner[i] >= 0)
                {
                  hbuff[i] = true;
                  c = hanging_corner[i];
                  num_parents = corner_num_hanging[i ^ c];
                  base_corner = corner_to_hanging[i ^ c];
                  for (j = 0; j < num_parents; ++j)
                    pbuff[j + 4 * i] = base_corner[j] ^ c;
                }
        }
      // Ghost elements.
      else
        {
          p4est_locidx_t idx = this->qtq -
            the_tmesh->num_local_octants ();
          
          for (i = 0; i < 8; ++i)
            {
              tbuff[i] = the_tmesh->ghost_data[40*idx + i];
              
              pbuff[4*i]     = the_tmesh->ghost_data[40*idx + 8  + 4*i];
              pbuff[4*i + 1] = the_tmesh->ghost_data[40*idx + 9  + 4*i];
              pbuff[4*i + 2] = the_tmesh->ghost_data[40*idx + 10 + 4*i];
              pbuff[4*i + 3] = the_tmesh->ghost_data[40*idx + 11 + 4*i];
              
              if (pbuff[4*i] != -1
                  || pbuff[4*i+1] != -1
                  || pbuff[4*i+2] != -1
                  || pbuff[4*i+3] != -1)
                hbuff[i] = true;
              else
                hbuff[i] = false;
            }
        }
    }
};

int
tmesh_3d::octant_t::parent (tmesh_3d::idx_t ip, tmesh_3d::idx_t in)
{
  assert (pbuff[ip + in * 4] >= 0);
  return tbuff[pbuff[ip + in * 4]];
};

int
tmesh_3d::octant_t::gparent (tmesh_3d::idx_t ip, tmesh_3d::idx_t in)
{  
  if (! is_ghost)
    {
      assert (pbuff[ip + in * 4] >= 0);

      return p8est_lnodes_global_index
        (the_tmesh->lnodes,
         static_cast<p4est_locidx_t>
         (tbuff[pbuff[ip + in * 4]]));
    }
  else
      return pbuff[4*in + ip];
};

tmesh_3d::idx_t
tmesh_3d::octant_t::e (idx_t i)
{  
  assert (i < 8);
  idx_t retval = NOT_ON_BOUNDARY;
  p8est_quadrant_t node;
  p8est_quadrant_corner_node (this->the_octant, i, &node);
  
  if (node.z == 0)
    retval = 4;
  else if (node.z == P8EST_ROOT_LEN)
    retval = 5;
  else if (node.y == 0)
    retval = 2;
  else if (node.y == P8EST_ROOT_LEN)
    retval = 3;
  else if (node.x == 0)
    retval = 0;
  else if (node.x == P8EST_ROOT_LEN)
    retval = 1;

  return retval;
};

tmesh_3d::neighbor_iterator
tmesh_3d::octant_t::begin_neighbor_sweep ()
{
  neighbor_iterator ni;
  
  if (this->the_tmesh->mesh == nullptr)
    this->the_tmesh->update ();
  
  p8est_mesh_face_neighbor_init (ni.face_neighbor,
                                 this->the_tmesh->p8est,
                                 this->the_tmesh->ghost,
                                 this->the_tmesh->mesh,
                                 this->get_tree_idx (),
                                 this->the_octant);
  
  p4est_topidx_t which_tree;
  p4est_locidx_t which_quad;
  int nface, nrank;
  
  p8est_quadrant_t * neighbor =
    p8est_mesh_face_neighbor_next (ni.face_neighbor, &which_tree,
                                   &which_quad, &nface, &nrank);
  
  ni.data = new octant_t (this->the_tmesh, which_tree, neighbor);
  
  p8est_tree_t *tree =
    p8est_tree_array_index (this->the_tmesh->p8est->trees,
                            which_tree);
  
  // If non-ghost.
  if (ni.face_neighbor->current_qtq <
      the_tmesh->num_local_octants ())
    {
      ni.data->forest_oct_idx = tree->quadrants_offset + which_quad;
      ni.data->tree_oct_idx = which_quad;
    }
  // If ghost.
  else
    {
      ni.data->is_ghost = true;
      ni.data->qtq = ni.face_neighbor->current_qtq;
      
      ni.data->forest_oct_idx =
        neighbor->p.piggy3.local_num +
        (the_tmesh->p8est->global_first_quadrant[nrank] -
         the_tmesh->p8est->global_first_quadrant[this->the_tmesh->rank]);
      
      ni.data->tree_oct_idx =
        ni.data->forest_oct_idx - tree->quadrants_offset;
    }
  
  ni.data->update (which_tree, neighbor);
  
  ni.face_idx = nface;
  return ni;
}

tmesh_3d::idx_t
tmesh_3d::octant_t::t (tmesh_3d::idx_t i)
{ return tbuff[i]; };

tmesh_3d::idx_t
tmesh_3d::octant_t::gt (tmesh_3d::idx_t i)
{
  if (! is_ghost)
    {
      if (the_tmesh->mesh == nullptr)
        the_tmesh->update ();
      
      return p8est_lnodes_global_index
        (the_tmesh->lnodes,
         static_cast<p4est_locidx_t> (tbuff[i]));
    }
  else
    return tbuff[i];
};

bool
tmesh_3d::octant_t::is_hanging (tmesh_3d::idx_t i)
{ return hbuff[i]; };


tmesh_3d::~tmesh_3d ()
{
    if (! (this->p8est  == nullptr)) p8est_destroy (this->p8est);
    if (! (this->conn   == nullptr)) p8est_connectivity_destroy (this->conn);
    if (! (this->lnodes == nullptr)) p8est_lnodes_destroy (this->lnodes);
    if (! (this->mesh   == nullptr)) p8est_mesh_destroy   (this->mesh);
    if (! (this->ghost  == nullptr)) p8est_ghost_destroy  (this->ghost);
    
    if (! (this->mirror_data == nullptr)) delete[] this->mirror_data;
    if (! (this->ghost_data  == nullptr)) delete[] this->ghost_data;
};

template <class p_type, class p_type_count,
          class t_type, class t_type_count>
static void
arrays2connectivity (const p_type *p_matrix_start,
                     const p_type_count num_vertices,
                     const t_type *t_matrix_start,
                     const t_type_count num_trees,
                     p8est_connectivity_t **conn)
{

  *conn =
    p8est_connectivity_new (num_vertices, num_trees, 0, 0, 0, 0);

  const p_type *p_iter = p_matrix_start;
  double *v_iter = &((*conn)->vertices[0]);

  auto p_matrix_end = p_matrix_start + num_vertices * 3;
  while (p_iter < p_matrix_end)
      *(v_iter++) = *(p_iter++);

  const t_type *t_iter = t_matrix_start;
  p4est_topidx_t *tv_iter = &((*conn)->tree_to_vertex[0]);

  auto t_matrix_end = t_matrix_start + num_trees * 9;
  while (t_iter < t_matrix_end)
    {
      for (int n = 0; n < 8; ++n)
        *(tv_iter++) = *(t_iter++) - 1;
      ++t_iter;
    }

  for (t_type tree = 0; tree < (*conn)->num_trees; ++tree)
    for (int face = 0; face < 6; ++face)
      {
        (*conn)->tree_to_tree[6 * tree + face] = tree;
        (*conn)->tree_to_face[6 * tree + face] = face;
      }

  assert (p8est_connectivity_is_valid (*conn));
  p8est_connectivity_complete (*conn);
};

/* Read a 3d p8est connectivity from a compressed octave
 * binary file. The file should contain a struct
 * named "msh" with fields "p" and "t". The former
 * should be the list of vertex coordinates while the
 * latter the list of element vertices.
 *
 * The standard p8est ordering is assumed for the nodes
 *
 * 7                     8
 *  +---------------------+
 *  |\                    |\
 *  | \                   | \
 *  |  \                  |  \
 *  |   \                 |   \
 *  |   5+---------------------+6
 *  |    |                |    |
 *  +----|----------------+    |
 *  3\   |               4 \   |
 *    \  |                  \  |
 *     \ |                   \ |
 *      \|                    \|
 *       +---------------------+
 *       1                     2
 */

static void
octbingz2connectivity
(const char *filename, p8est_connectivity_t **conn)
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

  arrays2connectivity (p_matrix.fortran_vec (),
                       num_vertices,
                       t_matrix.fortran_vec (),
                       num_trees, conn);
};


void
tmesh_3d::read_connectivity (const char *filename, int source)
{
  if (rank == source)
    octbingz2connectivity (filename, &conn);
  
  conn = p8est_connectivity_bcast (conn, source, comm);
  p8est = p8est_new (comm, conn, 0, NULL, this);
};

void
tmesh_3d::read_connectivity (const double *p,
                             const p4est_topidx_t num_vertices,
                             const p4est_topidx_t *t,
                             const p4est_topidx_t num_trees,
                             int source)
{
  if (rank == source)
    arrays2connectivity (p, num_vertices,
                         t, num_trees, &conn);
  
  conn = p8est_connectivity_bcast (conn, source, comm);
  p8est = p8est_new (comm, conn, 0, NULL, this);
};

void
tmesh_3d::save (const char *filename)
{ p8est_save (filename, p8est, 0); };

void
tmesh_3d::load (const char *filename)
{ p8est = p8est_load (filename, comm, 0, 0, this, &conn); };

void
tmesh_3d::vtk_export (const char *filename)
{
  p8est_vtk_context_t *context =
    p8est_vtk_context_new (p8est, filename);
  assert (context != nullptr);
  p8est_vtk_context_set_scale (context, 1.0);
  p8est_vtk_context_set_continuous (context, 1);
  context = p8est_vtk_write_header (context);
  context =
    p8est_vtk_write_cell_dataf (context, 1, 1, 1, 0, 0, 0, context);
  assert (p8est_vtk_write_footer (context) == 0);
};

void
tmesh_3d::update ()
{
  ghost  = p8est_ghost_new  (p8est, P8EST_CONNECT_FULL);
  lnodes = p8est_lnodes_new (p8est, ghost, 1);
  mesh   = p8est_mesh_new   (p8est, ghost, P8EST_CONNECT_FULL);
  
  update_ghosts ();
}

void
tmesh_3d::update_ghosts ()
{
  // Send mirror data.
  constexpr p4est_locidx_t chunk_len = 40;
  constexpr size_t data_size = sizeof (p4est_gloidx_t);
  
  p4est_locidx_t mirror_data_len =
    ghost->mirror_proc_offsets[size] * chunk_len;
  mirror_data = new p4est_gloidx_t[mirror_data_len];
  
  int mirror_end = 0;
  int mirror_begin = 0;
  
  std::vector<MPI_Request> req_s;
  p4est_locidx_t start, end, n_mirror;

  int tag = 0;
  int send_size = 0;
  
  // Loop over ranks.
  for (int i = 0; i < size; ++i)
    {
      start    = ghost->mirror_proc_offsets[i];
      end      = ghost->mirror_proc_offsets[i+1];
      n_mirror = end - start;
      
      p8est_quadrant_t * q;
      
      // Loop over mirrors.
      for (p4est_locidx_t j = start; j < end; ++j)
        {
          q =
            p8est_quadrant_array_index (&ghost->mirrors,
                                        ghost->mirror_proc_mirrors[j]);
          
          p8est_tree_t * tree =
            p8est_tree_array_index (p8est->trees, q->p.which_tree);
      
          idx_t global_idx =
            p8est->global_first_quadrant[rank] +
            q->p.piggy3.local_num;
          
          octant_t current_mirror (this, q->p.which_tree, q);
          current_mirror.forest_oct_idx = q->p.piggy3.local_num;
          current_mirror.tree_oct_idx =
            current_mirror.forest_oct_idx - tree->quadrants_offset;
          current_mirror.update (q->p.which_tree, q);
          
          for (int node = 0; node < 8; ++node)
            mirror_data[mirror_end++] = current_mirror.gt (node);
          
          for (int node = 0; node < 8; ++node)
            if (current_mirror.is_hanging (node))
              {
                mirror_data[mirror_end++] =
                  current_mirror.gparent (0, node);
                mirror_data[mirror_end++] =
                  current_mirror.gparent (1, node);
                mirror_data[mirror_end++] =
                  current_mirror.gparent (2, node);
                mirror_data[mirror_end++] =
                  current_mirror.gparent (3, node);
              }
            else
              {
                mirror_data[mirror_end++] = -1;
                mirror_data[mirror_end++] = -1;
                mirror_data[mirror_end++] = -1;
                mirror_data[mirror_end++] = -1;
              }
        }
      
      if (n_mirror > 0)
        {
          MPI_Request req;
          tag = rank + size * i;
          send_size = chunk_len * data_size * n_mirror;
          MPI_Isend (&(mirror_data[mirror_begin]), send_size,
                     MPI_CHAR, i, tag, comm, &req);
          req_s.push_back (req);
          
          /*std::cout << "Rank " << rank
                    << " is sending mirrors to rank "
                    << i << "." << std::endl;*/
        }
      
      mirror_begin = mirror_end;
    }
  
  // Receive ghost data.
  p4est_locidx_t ghosts_data_len =
    ghost->ghosts.elem_count * chunk_len;
  
  ghost_data = new p4est_gloidx_t[ghosts_data_len];
  int ghost_begin = 0;  
  p4est_locidx_t n_ghosts = 0;
  int recv_size = 0;
  // Loop over ranks.
  for (int i = 0; i < size; ++i)
    {
      start    = ghost->proc_offsets[i];
      end      = ghost->proc_offsets[i+1];
      n_ghosts = end - start;
      
      if (n_ghosts > 0)
        {
          MPI_Request req;
          tag = i + size * rank;
          recv_size = chunk_len * data_size * n_ghosts;
          MPI_Irecv (&(ghost_data[ghost_begin]), recv_size,
                     MPI_CHAR, i, tag, comm, &req);
          req_s.push_back (req);
          
          /*std::cout << "Rank " << rank
                    << " is receiving ghosts from rank "
                    << i << "." << std::endl;*/
        }
      
      ghost_begin += chunk_len * n_ghosts;
    }
  
  std::vector<MPI_Status> stats (req_s.size ());
  MPI_Waitall (req_s.size (), &(req_s[0]), &(stats[0]));
};

std::vector<int>
tmesh_3d::userint_replace (std::vector<int> old_userint)
{
  std::vector<int> new_userint;
  
  // Refinement.
  if (old_userint.size () == 1)
    {
      new_userint.resize (8);
      
      for (size_t i = 0; i < new_userint.size (); ++i)
        new_userint[i] = old_userint[0] - 1;
    }
  // Coarsening.
  else if (old_userint.size () == 8)
    {
      new_userint.resize (1);
      
      new_userint[0] = *std::min_element (old_userint.begin (),
                                          old_userint.end ()) + 1;
    }
  
  return new_userint;
}
