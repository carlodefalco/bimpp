/*! \file tmesh_3d.cpp
  \brief Interface for p8est library
*/

#include <array>
#include <iostream>
#include <octave_file_io.h>


#include <tmesh_3d.h>



double
tmesh_3d::quadrant_t::p (tmesh_3d::idx_t ii, tmesh_3d::idx_t jj)
{
  double retval = vxyz[3 * jj + ii];
  return (retval);
};

double
tmesh_3d::quadrant_t::centroid (tmesh_3d::idx_t ii)
{
  double retval = 0;

  for (tmesh_3d::idx_t c = 0; c < P8EST_CHILDREN; ++c)
    retval += this->p (ii, c);
  retval /= P8EST_CHILDREN;

  return (retval);
};

double
tmesh_3d::quadrant_t::face_centroid (tmesh_3d::idx_t ii, int jj)
{
  double retval = 0;

  for (tmesh_3d::idx_t c = 0; c < P8EST_HALF; ++c)
    retval += this->p (ii, p8est_face_corners[jj][c]);
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
static const int  zero = 0; /**< Constant zero. */
static const int  ones = P8EST_CHILDREN - 1; /**< One bit per dimension. */
static const int *corner_to_hanging[P8EST_CHILDREN];
static const int  corner_num_hanging[P8EST_CHILDREN] = {1, 2, 2, 4,
                                                        2, 4, 4, 1};

static int
lnodes_decode2 (p8est_lnodes_code_t face_code,
                int hanging_corner[P8EST_CHILDREN])
{
  if (face_code)
    {
      const int           c = (int) (face_code & ones);
      int                 i, h;
      int                 work = (int) (face_code >> P8EST_DIM);

      /* These two corners are never hanging by construction. */
      hanging_corner[c] = hanging_corner[c ^ ones] = -1;
      for (i = 0; i < P8EST_DIM; ++i)
        {
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
tmesh_3d::quadrant_iterator::operator++ ()
{

  p8est_t *p8 = data->the_tmesh->p8est;

  data->forest_quad_idx++;
  data->tree_quad_idx++;

  if (data->tree_quad_idx >= data->num_quadrants)
    {
      data->tree_idx++;

      if ((data->tree_idx) > (p8->last_local_tree))
        {
          this->data = nullptr;
          return;
        }
      data->tree_quad_idx = 0;

      data->tree = p8est_tree_array_index (p8->trees, data->tree_idx);
      data->tquadrants = &(data->tree)->quadrants;

      data->num_quadrants =
        (p4est_locidx_t) data->tquadrants->elem_count;
    }

  auto tmp = p8est_quadrant_array_index (data->tquadrants,
                                         data->tree_quad_idx);
  data->update (data->tree_idx, tmp);

};

void
tmesh_3d::neighbor_iterator::operator++ ()
{
  p4est_topidx_t which_tree;
  p4est_locidx_t which_quad;
  int nface, nrank;
  tmesh_3d *tmsh = data->the_tmesh;
  p8est_t *p8 = tmsh->p8est;

  p8est_quadrant_t * neighbor =
    p8est_mesh_face_neighbor_next (face_neighbor, &which_tree,
                                   &which_quad, &nface, &nrank);

  if (neighbor != nullptr)
    {
      p8est_tree_t * tree =
        p8est_tree_array_index (p8->trees,
                                which_tree);

      // If non-ghost.
      if (face_neighbor->current_qtq <
          tmsh->num_local_quadrants ())
        {
          data->is_ghost = false;
          data->qtq = -1;

          data->forest_quad_idx = tree->quadrants_offset + which_quad;
          data->tree_quad_idx = which_quad;
        }
      // If ghost.
      else
        {
          data->is_ghost = true;
          data->qtq = face_neighbor->current_qtq;

          data->forest_quad_idx =
            neighbor->p.piggy3.local_num +
            (p8->global_first_quadrant[nrank] -
             p8->global_first_quadrant[tmsh->rank]);

          data->tree_quad_idx = data->forest_quad_idx -
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
tmesh_3d::quadrant_t::update (p4est_topidx_t tree,
                              p8est_quadrant_t *q)
{
  p8est_quadrant_t node, parent;
  idx_t i, j;
  int hanging_corner[8];
  p8est_lnodes_t *ln = the_tmesh->lnodes;
  p4est_locidx_t lni;

  int c, h, n_parents;
  const int *base_corner;

  corner_to_hanging[0]        = &zero;
  corner_to_hanging[1]        = p8est_edge_corners[0];
  corner_to_hanging[2]        = p8est_edge_corners[4];
  corner_to_hanging[3]        = p8est_face_corners[4];
  corner_to_hanging[4]        = p8est_edge_corners[8];
  corner_to_hanging[ones - 2] = p8est_face_corners[2];
  corner_to_hanging[ones - 1] = p8est_face_corners[0];
  corner_to_hanging[ones]     = &ones;

  this->tree_idx     = tree;
  this->the_quadrant = q;

  for (i = 0; i < 8; ++i)
    {
      p8est_quadrant_corner_node (this->the_quadrant, i, &node);
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
              tbuff[i] = ln->element_nodes[8 * forest_quad_idx + i];
              hbuff[i] = 0;
              for (j = 0; j < 4; ++j)
                pbuff[4 * i + j] = -1;
            }

          bool any_hanging =
            lnodes_decode2 (ln->face_code[forest_quad_idx],
                            hanging_corner);
          if (any_hanging)
            for (i = 0; i < 8; ++i)
              if (hanging_corner[i] >= 0)
                {
                  c = hanging_corner[i];
                  n_parents = corner_num_hanging[i ^ c];
                  base_corner = corner_to_hanging[i ^ c];
                  hbuff[i] = n_parents;
                  for (j = 0; j < n_parents; ++j)
                    pbuff[4 * i + j] = base_corner[j] ^ c;
                }
        }
      // Ghost elements.
      else
        {
          p4est_locidx_t idx = this->qtq -
            the_tmesh->num_local_quadrants ();

          for (i = 0; i < 8; ++i)
            {
              tbuff[i] = the_tmesh->ghost_data[40 * idx + i];

              pbuff[4 * i]     =
                the_tmesh->ghost_data[40 * idx + 8  + 4 * i];
              pbuff[4 * i + 1] =
                the_tmesh->ghost_data[40 * idx + 9  + 4 * i];
              pbuff[4 * i + 2] =
                the_tmesh->ghost_data[40 * idx + 10 + 4 * i];
              pbuff[4 * i + 3] =
                the_tmesh->ghost_data[40 * idx + 11 + 4 * i];

              for (j = 0; j < 4; ++j)
                if (pbuff[4 * i + j] == -1)
                  break;
              hbuff[i] = j;
            }
        }
    }
};

int
tmesh_3d::quadrant_t::parent (tmesh_3d::idx_t ip, tmesh_3d::idx_t in)
{
  assert (pbuff[4 * in + ip] >= 0);
  return tbuff[pbuff[4 * in + ip]];
};

int
tmesh_3d::quadrant_t::gparent (tmesh_3d::idx_t ip, tmesh_3d::idx_t in)
{
  if (! is_ghost)
    {
      assert (pbuff[4 * in + ip] >= 0);

      return p8est_lnodes_global_index
        (the_tmesh->lnodes,
         static_cast<p4est_locidx_t>
         (tbuff[pbuff[4 * in + ip]]));
    }
  else
    return pbuff[4 * in + ip];
};

tmesh_3d::idx_t
tmesh_3d::quadrant_t::e (idx_t i)
{
  assert (i < 8);
  idx_t retval = NOT_ON_BOUNDARY;
  p8est_quadrant_t node;
  p8est_quadrant_corner_node (this->the_quadrant, i, &node);

  if      (node.z == 0)
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
tmesh_3d::quadrant_t::begin_neighbor_sweep ()
{
  neighbor_iterator ni;

  if (this->the_tmesh->mesh == nullptr)
    this->the_tmesh->update ();

  p8est_mesh_face_neighbor_init (ni.face_neighbor,
                                 this->the_tmesh->p8est,
                                 this->the_tmesh->ghost,
                                 this->the_tmesh->mesh,
                                 this->get_tree_idx (),
                                 this->the_quadrant);

  p4est_topidx_t which_tree;
  p4est_locidx_t which_quad;
  int nface, nrank;

  p8est_quadrant_t * neighbor =
    p8est_mesh_face_neighbor_next (ni.face_neighbor, &which_tree,
                                   &which_quad, &nface, &nrank);

  ni.data = new quadrant_t (this->the_tmesh, which_tree, neighbor);

  p8est_tree_t *tree =
    p8est_tree_array_index (this->the_tmesh->p8est->trees,
                            which_tree);

  // If non-ghost.
  if (ni.face_neighbor->current_qtq <
      the_tmesh->num_local_quadrants ())
    {
      ni.data->forest_quad_idx = tree->quadrants_offset + which_quad;
      ni.data->tree_quad_idx = which_quad;
    }
  // If ghost.
  else
    {
      ni.data->is_ghost = true;
      ni.data->qtq = ni.face_neighbor->current_qtq;

      ni.data->forest_quad_idx =
        neighbor->p.piggy3.local_num +
        (the_tmesh->p8est->global_first_quadrant[nrank] -
         the_tmesh->p8est->global_first_quadrant[this->the_tmesh->rank]);

      ni.data->tree_quad_idx =
        ni.data->forest_quad_idx - tree->quadrants_offset;
    }

  ni.data->update (which_tree, neighbor);

  ni.face_idx = nface;
  return ni;
}

tmesh_3d::idx_t
tmesh_3d::quadrant_t::t (tmesh_3d::idx_t i)
{ return tbuff[i]; };

tmesh_3d::idx_t
tmesh_3d::quadrant_t::gt (tmesh_3d::idx_t i)
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
tmesh_3d::quadrant_t::is_hanging (tmesh_3d::idx_t i)
{ return hbuff[i] > 0; };

int
tmesh_3d::quadrant_t::num_parents (tmesh_3d::idx_t i)
{ return hbuff[i]; };

tmesh_3d::~tmesh_3d ()
{
  if (! (this->p8est       == nullptr))
    p8est_destroy (this->p8est);
  if (! (this->conn        == nullptr))
    p8est_connectivity_destroy (this->conn);
  if (! (this->lnodes      == nullptr))
    p8est_lnodes_destroy (this->lnodes);
  if (! (this->mesh        == nullptr))
    p8est_mesh_destroy   (this->mesh);
  if (! (this->ghost       == nullptr))
    p8est_ghost_destroy  (this->ghost);
  if (! (this->mirror_data == nullptr))
    delete[] this->mirror_data;
  if (! (this->ghost_data  == nullptr))
    delete[] this->ghost_data;
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
 * The standard p8est ordering is assumed for the nodes (1 based !!!)
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

  int flag_open = octave_io_open (filename, m, &m);
  assert (flag_open == 0);

  int flag_load = octave_load ("msh", tmp);
  assert (flag_load == 0);
  
  Matrix p_matrix =
    tmp.scalar_map_value ().contents ("p").matrix_value ();

  Array<octave_idx_type> t_matrix =
    tmp.scalar_map_value ().contents ("t").octave_idx_type_vector_value (false, true, true);

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
  p8est = p8est_new (comm, conn, sizeof (tmesh_3d::data_t), nullptr, this);
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
  p8est = p8est_new (comm, conn, sizeof (tmesh_3d::data_t), nullptr, this);
};

void
tmesh_3d::save (const char *filename)
{ p8est_save (filename, p8est, 0); };

void
tmesh_3d::load (const char *filename)
{ p8est = p8est_load (filename, comm, sizeof (tmesh_3d::data_t), 0, 
                      this, &conn); };

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

  int flag = p8est_vtk_write_footer (context);
  assert (flag == 0);
};

template<class T>
void
octbin_export_tmpl (tmesh_3d *THIS, const char* basename, const T& f,
                    const ordering& ord)
{
  //  assert (f.size () == num_global_nodes ());

  std::vector<double> p (3 * THIS->num_owned_nodes ());
  std::vector<double> f_loc (THIS->num_owned_nodes ());

  Array<octave_idx_type>
    oct_t (dim_vector (8, THIS->num_local_quadrants ()), 0);
  octave_idx_type *t = oct_t.fortran_vec ();

  octave_idx_type ij = 0;
  for (auto quadrant = THIS->begin_quadrant_sweep ();
       quadrant != THIS->end_quadrant_sweep ();
       ++quadrant)
    {
      ij = 0;
      for (int ii = 0; ii < 8; ++ii)
        {
          if (! quadrant->is_hanging (ii))
            if (quadrant->t (ii) < THIS->num_owned_nodes ())
	            {
	              for (int jj = 0; jj < 3; ++jj)
                  p[3 * quadrant->t (ii) + jj] = quadrant->p (jj, ii);
                f_loc[quadrant->t (ii)] = f[ord(quadrant->gt (ii))];
                t[8 * quadrant->get_forest_quad_idx () + (ij++)] =
                  quadrant->t (ii);
	            }
          	else
	            {
	              for (int jj = 0; jj < 3; ++jj)
                  p.push_back (quadrant->p (jj, ii));
                f_loc.push_back (f[ord(quadrant->gt (ii))]);
                t[8 * quadrant->get_forest_quad_idx () + (ij++)] =
                  f_loc.size () - 1;
	            }
          else
            {
              for (int jj = 0; jj < 3; ++jj)
                p.push_back (quadrant->p (jj, ii));
              int pp;
              double fbuff = 0;
              for (pp = 0; pp < quadrant->num_parents (ii); ++pp)
                fbuff += f [ord(quadrant->gparent (pp, ii))];
              f_loc.push_back (fbuff / pp);
              t[8 * quadrant->get_forest_quad_idx () + (ij++)] =
                f_loc.size () - 1;
            }
        }
    }
 
  Matrix oct_p (3, p.size () / 3, 0.0);
  ColumnVector oct_f (f_loc.size (), 0.0);

  std::copy_n (p.begin (), p.size (), oct_p.fortran_vec ());
  std::copy_n (f_loc.begin (), f_loc.size (), oct_f.fortran_vec ());

  octave_scalar_map the_map;
  the_map.assign ("p", oct_p);
  the_map.assign ("f", oct_f);
  the_map.assign ("t", oct_t);

  octave_io_mode m = gz_write_mode;

  // Define filename.
  char filename[255] = "";
  sprintf (filename, "%s_%4.4d.octbin.gz", basename, THIS->rank);

  // Save to filename.
  int flag_open = octave_io_open (filename, m, &m);
  assert (flag_open == 0);
  
  int flag_save = octave_save ("msh", octave_value (the_map));
  assert (flag_save == 0);

  int flag_close = octave_io_close ();
  assert (flag_close == 0);
  
};

void
tmesh_3d::octbin_export (const char * filename,
               const distributed_vector & f, 
               const ordering& ord)
{
	octbin_export_tmpl (this, filename, f, ord);
};

void
tmesh_3d::octbin_export (const char * filename,
                         const std::vector<double> & f,
                         const ordering& ord)
{
	octbin_export_tmpl (this, filename, f, ord);
}	

void
tmesh_3d::octbin_export_quadrant (const char * basename,
                                  const std::vector<double> & f)
{
  assert (f.size () == num_local_quadrants ());

  ColumnVector oct_f (f.size (), 0.0);

  std::copy_n (f.begin (), f.size (), oct_f.fortran_vec ());

  octave_scalar_map the_map;
  the_map.assign ("f", oct_f);
  
  octave_io_mode m = gz_write_mode;

  // Define filename.
  char filename[255] = "";
  sprintf (filename, "%s_%4.4d.octbin.gz", basename, rank);

  // Save to filename.
  int flag_open = octave_io_open (filename, m, &m);
  assert (flag_open == 0);
  
  int flag_save = octave_save ("msh", octave_value (the_map));
  assert (flag_save == 0);

  int flag_close = octave_io_close ();
  assert (flag_close == 0);
};

void
tmesh_3d::quadrant_iterator::reset ()
{
  if (data != nullptr)
    {
      p8est_t *p8           = data->the_tmesh->p8est;
      data->tree_idx        = p8->first_local_tree;
      data->tree_quad_idx   = 0;
      data->forest_quad_idx = 0;

      data->is_ghost = false;
      data->qtq = -1;

      if (data->tree_idx != -1)
        {
          data->tree          =
            p8est_tree_array_index (p8->trees, data->tree_idx);
          data->tquadrants    = &(data->tree->quadrants);
          data->num_quadrants =
            (p4est_locidx_t) data->tquadrants->elem_count;
        }
      else
        {
          data->tree          = nullptr;
          data->tquadrants    = nullptr;
          data->num_quadrants = 0;
        }

      if (data->num_quadrants > 0)
        {
          auto tmp =
            p8est_quadrant_array_index (data->tquadrants,
                                        data->forest_quad_idx);
          data->update (data->tree_idx, tmp);
        }
      else
        data = nullptr;
    }
};

tmesh_3d::quadrant_iterator
tmesh_3d::begin_quadrant_sweep ()
{
  if (this->mesh == nullptr)
    this->update ();

  quadrant_iterator qi (&current_quadrant);
  qi.reset ();
  return qi;
};

void
tmesh_3d::set_metrics_marker
(std::function<double (tmesh_3d::quadrant_iterator)> estimator,
 double tol, int max_depth, int n_refine, int n_coarsen)
{
  this->metrics_max_depth = max_depth;

  tmesh_3d::data_t * data;
  int hxhat_hx = 0;

  for (auto quadrant = this->begin_quadrant_sweep ();
       quadrant != this->end_quadrant_sweep (); ++quadrant)
    {
      set_interpolation_matrix (quadrant);

      hxhat_hx = static_cast<int> (std::round(std::log2 (estimator (quadrant)
                          * std::sqrt (this->num_global_quadrants ()) / tol)));

      if (hxhat_hx >= 0)
        hxhat_hx = std::max(0, hxhat_hx - n_refine);
      else
        hxhat_hx = std::min(0, hxhat_hx - n_coarsen);

      data = 
        static_cast<tmesh_3d::data_t *> (quadrant->the_quadrant->p.user_data);
      
      data->refine_count = 
        std::min (std::max (-max_depth, hxhat_hx), max_depth);
    }

  return;
}

void
tmesh_3d::refine (int recursive, int partforcoarsen, int balance)
{
  quadrant_iterator qi (&current_quadrant);
  qi.reset ();

  p8est_refine_ext (p8est, recursive, -1, refine_callback,
                    nullptr, replace_callback);

  if (balance)
    p8est_balance_ext (p8est, P8EST_CONNECT_EDGE, nullptr, replace_callback);

  p8est_partition (p8est, partforcoarsen, nullptr);

  if (! (lnodes == nullptr)) p8est_lnodes_destroy (lnodes);
  lnodes = nullptr;

  if (! (mesh == nullptr)) p8est_mesh_destroy (mesh);
  mesh = nullptr;

  if (! (ghost == nullptr)) p8est_ghost_destroy (ghost);
  ghost = nullptr;
}

void
tmesh_3d::metrics_refine (idx_t max_elems)
{
  int recursive = 0;
  int partforcoarsen = 1;
  int balance = 0;

  for (int i = 0; i < metrics_max_depth; ++i)
    {
      coarsen (recursive, partforcoarsen, balance);

      // Prevent large meshes.
      if (max_elems <= 0 || this->num_global_quadrants () < max_elems)
        refine (recursive, partforcoarsen, balance);
    }

  p8est_balance_ext (p8est, P8EST_CONNECT_EDGE, nullptr, replace_callback);
}

void
tmesh_3d::coarsen (int recursive, int partforcoarsen, int balance)
{
  p8est_coarsen_ext (p8est, recursive, 0, coarsen_callback,
                      nullptr, replace_callback);

  if (balance)
    p8est_balance_ext (p8est, P8EST_CONNECT_EDGE, nullptr, replace_callback);

  p8est_partition (p8est, partforcoarsen, nullptr);

  if (! (lnodes == nullptr)) p8est_lnodes_destroy (lnodes);
  lnodes = nullptr;

  if (! (mesh == nullptr)) p8est_mesh_destroy (mesh);
  mesh = nullptr;

  if (! (ghost == nullptr)) p8est_ghost_destroy (ghost);
  ghost = nullptr;
};

void
tmesh_3d::update ()
{
  ghost  = p8est_ghost_new  (p8est, P8EST_CONNECT_FULL);
  lnodes = p8est_lnodes_new (p8est, ghost, 1);
  mesh   = p8est_mesh_new   (p8est, ghost, P8EST_CONNECT_EDGE);

  update_ghosts ();
}

void
tmesh_3d::update_ghosts ()
{
  // deallocate previously built storage arrays
  if (! (this->mirror_data == nullptr))
    delete[] this->mirror_data;
  if (! (this->ghost_data  == nullptr))
    delete[] this->ghost_data;
  
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
      end      = ghost->mirror_proc_offsets[i + 1];
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

          quadrant_t current_mirror (this, q->p.which_tree, q);
          current_mirror.forest_quad_idx = q->p.piggy3.local_num;
          current_mirror.tree_quad_idx =
            current_mirror.forest_quad_idx - tree->quadrants_offset;
          current_mirror.update (q->p.which_tree, q);

          for (int node = 0; node < 8; ++node)
            mirror_data[mirror_end++] = current_mirror.gt (node);

          for (int node = 0; node < 8; ++node)
            for (int pp = 0; pp < 4; ++pp)
              mirror_data[mirror_end++] =
                current_mirror.num_parents (node) > pp ?
                current_mirror.gparent (pp, node) : -1;
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
      end      = ghost->proc_offsets[i + 1];
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

std::vector<tmesh_3d::data_t>
tmesh_3d::user_data_replace (std::vector<tmesh_3d::data_t *> old_user_data)
{
  std::vector<tmesh_3d::data_t> new_user_data;

  // Refinement.
  if (old_user_data.size () == 1)
    {
      new_user_data.resize (8);

      for (size_t i = 0; i < new_user_data.size (); ++i)
        {
          // Decrease refine_count.
          new_user_data[i].refine_count =
            old_user_data[0]->refine_count - 1;

          // Determine interpolation indices.
          new_user_data[i].interp_idx =
            old_user_data[0]->interp_idx;
          
          // Compute local interpolation matrix.
          std::array<std::array<double, 8>, 8> loc_interp;
          
          if (i == 0)
            loc_interp =
              {
                1,     0,     0,     0,     0,     0,     0,     0,
                0.5,   0.5,   0,     0,     0,     0,     0,     0,
                0.5,   0,     0.5,   0,     0,     0,     0,     0,
                0.25,  0.25,  0.25,  0.25,  0,     0,     0,     0,
                0.5,   0,     0,     0,     0.5,   0,     0,     0,
                0.25,  0.25,  0,     0,     0.25,  0.25,  0,     0,
                0.25,  0,     0.25,  0,     0.25,  0,     0.25,  0,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125
              };
          else if (i == 1)
            loc_interp =
              {
                0.5,   0.5,   0,     0,     0,     0,     0,     0,
                0,     1,     0,     0,     0,     0,     0,     0,
                0.25,  0.25,  0.25,  0.25,  0,     0,     0,     0,
                0,     0.5,   0,     0.5,   0,     0,     0,     0,
                0.25,  0.25,  0,     0,     0.25,  0.25,  0,     0,
                0,     0.5,   0,     0,     0,     0.5,   0,     0,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                0,     0.25,  0,     0.25,  0,     0.25,  0,     0.25
              };
          else if (i == 2)
            loc_interp =
              {
                0.5,   0,     0.5,   0,     0,     0,     0,     0,
                0.25,  0.25,  0.25,  0.25,  0,     0,     0,     0,
                0,     0,     1,     0,     0,     0,     0,     0,
                0,     0,     0.5,   0.5,   0,     0,     0,     0,
                0.25,  0,     0.25,  0,     0.25,  0,     0.25,  0,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                0,     0,     0.5,   0,     0,     0,     0.5,   0,
                0,     0,     0.25,  0.25,  0,     0,     0.25,  0.25
              };
          else if (i == 3)
            loc_interp =
              {
                0.25,  0.25,  0.25,  0.25,  0,     0,     0,     0,
                0,     0.5,   0,     0.5,   0,     0,     0,     0,
                0,     0,     0.5,   0.5,   0,     0,     0,     0,
                0,     0,     0,     1,     0,     0,     0,     0,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                0,     0.25,  0,     0.25,  0,     0.25,  0,     0.25,
                0,     0,     0.25,  0.25,  0,     0,     0.25,  0.25,
                0,     0,     0,     0.5,   0,     0,     0,     0.5
              };
          else if (i == 4)
            loc_interp =
              {
                0.5,   0,     0,     0,     0.5,   0,     0,     0,
                0.25,  0.25,  0,     0,     0.25,  0.25,  0,     0,
                0.25,  0,     0.25,  0,     0.25,  0,     0.25,  0,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                0,     0,     0,     0,     1,     0,     0,     0,
                0,     0,     0,     0,     0.5,   0.5,   0,     0,
                0,     0,     0,     0,     0.5,   0,     0.5,   0,
                0,     0,     0,     0,     0.25,  0.25,  0.25,  0.25
              };
          else if (i == 5)
            loc_interp =
              {
                0.25,  0.25,  0,     0,     0.25,  0.25,  0,     0,
                0,     0.5,   0,     0,     0,     0.5,   0,     0,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                0,     0.25,  0,     0.25,  0,     0.25,  0,     0.25,
                0,     0,     0,     0,     0.5,   0.5,   0,     0,
                0,     0,     0,     0,     0,     1,     0,     0,
                0,     0,     0,     0,     0.25,  0.25,  0.25,  0.25,
                0,     0,     0,     0,     0,     0.5,   0,     0.5
              };
          else if (i == 6)
            loc_interp =
              {
                0.25,  0,     0.25,  0,     0.25,  0,     0.25,  0,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                0,     0,     0.5,   0,     0,     0,     0.5,   0,
                0,     0,     0.25,  0.25,  0,     0,     0.25,  0.25,
                0,     0,     0,     0,     0.5,   0,     0.5,   0,
                0,     0,     0,     0,     0.25,  0.25,  0.25,  0.25,
                0,     0,     0,     0,     0,     0,     1,     0,
                0,     0,     0,     0,     0,     0,     0.5,   0.5
              };
          else if (i == 7)
            loc_interp =
              {
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                0,     0.25,  0,     0.25,  0,     0.25,  0,     0.25,
                0,     0,     0.25,  0.25,  0,     0,     0.25,  0.25,
                0,     0,     0,     0.5,   0,     0,     0,     0.5,
                0,     0,     0,     0,     0.25,  0.25,  0.25,  0.25,
                0,     0,     0,     0,     0,     0.5,   0,     0.5,                
                0,     0,     0,     0,     0,     0,     0.5,   0.5,
                0,     0,     0,     0,     0,     0,     0,     1
              };
          
          // Multiply by parent interpolation matrix.
          new_user_data[i].interp_coeff = {0};
          
          for (int row = 0; row < 8; ++row)
            for (int col = 0; col < 8; ++col)
              for (int k = 0; k < 8; ++k)
                new_user_data[i].interp_coeff[row][col] +=
                  loc_interp[row][k] * old_user_data[0]->interp_coeff[k][col];
        }
    }
  // Coarsening.
  else if (old_user_data.size () == 8)
    {
      new_user_data.resize (1);

      // Increase refine_count.
      std::array<int, 8> ref_counts =
        {
          old_user_data[0]->refine_count,
          old_user_data[1]->refine_count,
          old_user_data[2]->refine_count,
          old_user_data[3]->refine_count,
          old_user_data[4]->refine_count,
          old_user_data[5]->refine_count,
          old_user_data[6]->refine_count,
          old_user_data[7]->refine_count
        };
      
      new_user_data[0].refine_count =
        *std::max_element (ref_counts.begin (),
                           ref_counts.end ()) + 1;

      // Replace interpolation matrix.
      new_user_data[0].interp_coeff = {0};
      
      for (int row = 0; row < 8; ++row)
        {
          // If coarsening, then (due to balancing)
          // the parent indices have a "1" entry.
          for (int col = 0; col < 8; ++col)
            {
              if (old_user_data[row]->interp_coeff[row][col] == 1)
                {
                  new_user_data[0].interp_idx[row] = 
                    old_user_data[row]->interp_idx[col];
                  
                  // Insert diagonal entry.
                  new_user_data[0].interp_coeff[row][row] = 1;
                  break;
                }
            }
        }
    }
  
  return new_user_data;
}

int
tmesh_3d::refine_callback (p8est_t* p8, p4est_topidx_t tt,
                           p8est_quadrant_t* qq)
{
  tmesh_3d::data_t * data =
    static_cast<tmesh_3d::data_t *> (qq->p.user_data);

  return (data->refine_count > 0);
};

int
tmesh_3d::coarsen_callback (p8est_t* p8, p4est_topidx_t tt,
                            p8est_quadrant_t* qq [])
{
  return (static_cast<tmesh_3d::data_t *> (qq[0]->p.user_data)->refine_count < 0 
      && static_cast<tmesh_3d::data_t *> (qq[1]->p.user_data)->refine_count < 0
      && static_cast<tmesh_3d::data_t *> (qq[2]->p.user_data)->refine_count < 0
      && static_cast<tmesh_3d::data_t *> (qq[3]->p.user_data)->refine_count < 0
      && static_cast<tmesh_3d::data_t *> (qq[4]->p.user_data)->refine_count < 0
      && static_cast<tmesh_3d::data_t *> (qq[5]->p.user_data)->refine_count < 0
      && static_cast<tmesh_3d::data_t *> (qq[6]->p.user_data)->refine_count < 0
      && static_cast<tmesh_3d::data_t *> (qq[7]->p.user_data)->refine_count < 0);
};

void
tmesh_3d::replace_callback (p8est_t * p8,
                            p4est_topidx_t tt,
                            int num_outgoing,
                            p8est_quadrant_t * outgoing[],
                            int num_incoming,
                            p8est_quadrant_t * incoming[])
{
  tmesh_3d *tm = reinterpret_cast<tmesh_3d*> (p8->user_pointer);

  std::vector<tmesh_3d::data_t *> old_user_data (num_outgoing);

  for (size_t i = 0; i < num_outgoing; ++i)
    old_user_data[i] = 
      static_cast<tmesh_3d::data_t *> (outgoing[i]->p.user_data);

  std::vector<tmesh_3d::data_t> new_user_data =
    tm->replace_fun (old_user_data);

  for (size_t i = 0; i < num_incoming; ++i)
    *(static_cast<tmesh_3d::data_t *> (incoming[i]->p.user_data)) =
      new_user_data[i];

  return;
};

void
tmesh_3d::set_interpolation_matrix (tmesh_3d::quadrant_iterator & q)
{
  // Create interpolation map.
  std::map<idx_t,
           std::vector<std::pair<int, double>>> interp_map;
  
  for (int node = 0; node < 8; ++node)
    {
      if (! q->is_hanging (node))
        interp_map[q->gt (node)].push_back
          (std::make_pair(node, 1));
      else
        {
          int np = q->num_parents(node);
          for (int pp = 0; pp < np; ++pp)
            interp_map[q->gparent (pp, node)].push_back
              (std::make_pair(node, 1/np));
        }
    }
  
  // Copy interp_map into user_data.
  tmesh_3d::data_t * data =
    static_cast<tmesh_3d::data_t *> (q->the_quadrant->p.user_data);
  
  data->interp_idx = {0};
  data->interp_coeff = {0};

  int col = 0;
  for (auto map_el = interp_map.begin ();
       map_el != interp_map.end ();
       ++col, ++map_el)
    {
      data->interp_idx[col] =
        map_el->first;
      
      for (auto vec_entry : map_el->second)
        data->interp_coeff[vec_entry.first][col] = vec_entry.second;
    }
};
