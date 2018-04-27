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
  double retval = vxyz[3 * jj + ii];
  return (retval);
};

double
tmesh::quadrant_t::centroid (tmesh::idx_t ii)
{
  double retval = 0;

  if (ii == 0)
    retval = 0.5 * (this->p (0, 0) + this->p (0, 1));
  else if (ii == 1)
    retval = 0.5 * (this->p (1, 0) + this->p (1, 2));

  return (retval);
};

/** Decode the information from p4est_lnodes_t for a given element.
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
static const int   *corner_to_hanging[4];
static const int    corner_num_hanging[4] = { 1, 2, 2, 1 };

static int
lnodes_decode2 (p4est_lnodes_code_t face_code,
                int hanging_corner[P4EST_CHILDREN])
{
  if (face_code)
    {
      const int           c = (int) (face_code & ones);
      int                 i, h;
      int                 work = (int) (face_code >> P4EST_DIM);

      /* These two corners are never hanging by construction. */
      hanging_corner[c] = hanging_corner[c ^ ones] = -1;
      for (i = 0; i < P4EST_DIM; ++i)
        {
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
tmesh::quadrant_iterator::operator++ ()
{

  p4est_t *p4 = data->the_tmesh->p4est;

  data->forest_quad_idx++;
  data->tree_quad_idx++;

  if (data->tree_quad_idx >= data->num_quadrants)
    {
      data->tree_idx++;

      if ((data->tree_idx) > (p4->last_local_tree))
        {
          this->data = nullptr;
          return;
        }
      data->tree_quad_idx = 0;

      data->tree = p4est_tree_array_index (p4->trees, data->tree_idx);
      data->tquadrants = &(data->tree)->quadrants;

      data->num_quadrants =
        (p4est_locidx_t) data->tquadrants->elem_count;
    }

  auto tmp = p4est_quadrant_array_index (data->tquadrants,
                                         data->tree_quad_idx);
  data->update (data->tree_idx, tmp);

};

void
tmesh::neighbor_iterator::operator++ ()
{
  p4est_topidx_t which_tree;
  p4est_locidx_t which_quad;
  int nface, nrank;
  tmesh *tmsh = data->the_tmesh;
  p4est_t *p4 = tmsh->p4est;

  p4est_quadrant_t * neighbor =
    p4est_mesh_face_neighbor_next (face_neighbor, &which_tree,
                                   &which_quad, &nface, &nrank);

  if (neighbor != nullptr)
    {
      p4est_tree_t * tree =
        p4est_tree_array_index (p4->trees,
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
            (p4->global_first_quadrant[nrank] -
             p4->global_first_quadrant[tmsh->rank]);

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
tmesh::quadrant_t::update (p4est_topidx_t tree,
                           p4est_quadrant_t *q)
{
  p4est_quadrant_t node, parent;
  idx_t i, j;
  int hanging_corner[4];
  p4est_lnodes_t *ln = the_tmesh->lnodes;
  p4est_locidx_t lni;

  int c, h, num_parents;
  const int *base_corner;

  corner_to_hanging[0]        = &zero;
  corner_to_hanging[ones - 2] = p4est_face_corners[2];
  corner_to_hanging[ones - 1] = p4est_face_corners[0];
  corner_to_hanging[ones]     = &ones;

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
      // Non-ghost elements.
      if (! is_ghost)
        {
          for (i = 0; i < 4; ++i)
            {
              tbuff[i] = ln->element_nodes[4 * forest_quad_idx + i];
              hbuff[i] = false;
              pbuff[2 * i] = -1;
              pbuff[2 * i + 1] = -1;
            }

          bool any_hanging =
            lnodes_decode2 (ln->face_code[forest_quad_idx],
                            hanging_corner);
          if (any_hanging)
            for (i = 0; i < 4; ++i)
              if (hanging_corner[i] >= 0)
                {
                  hbuff[i] = true;
                  c = hanging_corner[i];
                  num_parents = corner_num_hanging[i ^ c];
                  base_corner = corner_to_hanging[i ^ c];
                  for (j = 0; j < num_parents; ++j)
                    pbuff[j + 2 * i] = base_corner[j] ^ c;
                }
        }
      // Ghost elements.
      else
        {
          p4est_locidx_t idx = this->qtq -
            the_tmesh->num_local_quadrants ();

          for (i = 0; i < 4; ++i)
            {
              tbuff[i] = the_tmesh->ghost_data[12 * idx + i];

              pbuff[2 * i]   = the_tmesh->ghost_data[12 * idx + 4 + 2 * i];
              pbuff[2 * i + 1] = the_tmesh->ghost_data[12 * idx + 5 + 2 * i];

              if (pbuff[2 * i] != -1
                  || pbuff[2 * i + 1] != -1)
                hbuff[i] = true;
              else
                hbuff[i] = false;
            }
        }
    }
};

int
tmesh::quadrant_t::parent (tmesh::idx_t ip, tmesh::idx_t in)
{
  assert (pbuff[ip + in * 2] >= 0);
  return tbuff[pbuff[ip + in * 2]];
};

int
tmesh::quadrant_t::gparent (tmesh::idx_t ip, tmesh::idx_t in)
{
  if (! is_ghost)
    {
      assert (pbuff[ip + in * 2] >= 0);

      return p4est_lnodes_global_index
        (the_tmesh->lnodes,
         static_cast<p4est_locidx_t>
         (tbuff[pbuff[ip + in * 2]]));
    }
  else
    return pbuff[2 * in + ip];
};

tmesh::idx_t
tmesh::quadrant_t::e (idx_t i)
{
  assert (i < 4);
  idx_t retval = NOT_ON_BOUNDARY;
  p4est_quadrant_t node;
  p4est_quadrant_corner_node (this->the_quadrant, i, &node);

  if (node.y == 0)
    retval = 2;
  else if (node.y == P4EST_ROOT_LEN)
    retval = 3;
  else if (node.x == 0)
    retval = 0;
  else if (node.x == P4EST_ROOT_LEN)
    retval = 1;

  return retval;
};

tmesh::neighbor_iterator
tmesh::quadrant_t::begin_neighbor_sweep ()
{
  neighbor_iterator ni;

  if (this->the_tmesh->mesh == nullptr)
    this->the_tmesh->update ();

  p4est_mesh_face_neighbor_init (ni.face_neighbor,
                                 this->the_tmesh->p4est,
                                 this->the_tmesh->ghost,
                                 this->the_tmesh->mesh,
                                 this->get_tree_idx (),
                                 this->the_quadrant);

  p4est_topidx_t which_tree;
  p4est_locidx_t which_quad;
  int nface, nrank;

  p4est_quadrant_t * neighbor =
    p4est_mesh_face_neighbor_next (ni.face_neighbor, &which_tree,
                                   &which_quad, &nface, &nrank);

  ni.data = new quadrant_t (this->the_tmesh, which_tree, neighbor);

  p4est_tree_t *tree =
    p4est_tree_array_index (this->the_tmesh->p4est->trees,
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
        (the_tmesh->p4est->global_first_quadrant[nrank] -
         the_tmesh->p4est->global_first_quadrant[this->the_tmesh->rank]);

      ni.data->tree_quad_idx =
        ni.data->forest_quad_idx - tree->quadrants_offset;
    }

  ni.data->update (which_tree, neighbor);

  ni.face_idx = nface;
  return ni;
}

tmesh::idx_t
tmesh::quadrant_t::t (tmesh::idx_t i)
{ return tbuff[i]; };

tmesh::idx_t
tmesh::quadrant_t::gt (tmesh::idx_t i)
{
  if (! is_ghost)
    {
      if (the_tmesh->mesh == nullptr)
        the_tmesh->update ();

      return p4est_lnodes_global_index
        (the_tmesh->lnodes,
         static_cast<p4est_locidx_t> (tbuff[i]));
    }
  else
    return tbuff[i];
};

bool
tmesh::quadrant_t::is_hanging (tmesh::idx_t i)
{ return hbuff[i]; };


tmesh::~tmesh ()
{
  if (! (this->p4est  == nullptr)) p4est_destroy (this->p4est);
  if (! (this->conn   == nullptr)) p4est_connectivity_destroy (
                                                               this->conn);
  if (! (this->lnodes == nullptr)) p4est_lnodes_destroy (this->lnodes);
  if (! (this->mesh   == nullptr)) p4est_mesh_destroy   (this->mesh);
  if (! (this->ghost  == nullptr)) p4est_ghost_destroy  (this->ghost);

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
                     p4est_connectivity_t **conn)
{

  std::array<t_type, 4> v = {0, 0, 0, 0};
  int face = 0;

  *conn =
    p4est_connectivity_new (num_vertices, num_trees, 0, 0);

  const p_type *p_iter = p_matrix_start;
  double *v_iter = &((*conn)->vertices[0]);

  auto p_matrix_end = p_matrix_start + num_vertices * 2;
  while (p_iter < p_matrix_end)
    {
      *(v_iter++) = *(p_iter++);
      *(v_iter++) = *(p_iter++);
      *(v_iter++) = 0;
    }

  const t_type *t_iter = t_matrix_start;
  p4est_topidx_t *tv_iter = &((*conn)->tree_to_vertex[0]);

  auto t_matrix_end = t_matrix_start + num_trees * 5;
  while (t_iter < t_matrix_end)
    {

      for (int n = 0; n < 4; ++n)
        v[n] = *(t_iter++);
      ++t_iter;

      *(tv_iter++) = v[0] - static_cast<t_type> (1);
      *(tv_iter++) = v[1] - static_cast<t_type> (1);
      *(tv_iter++) = v[3] - static_cast<t_type> (1);
      *(tv_iter++) = v[2] - static_cast<t_type> (1);

    }

  for (t_type_count tree = 0; tree < num_trees; ++tree)
    for (face = 0; face < 4; ++face)
      {
        (*conn)->tree_to_tree[4 * tree + face] = tree;
        (*conn)->tree_to_face[4 * tree + face] = face;
      }

  assert (p4est_connectivity_is_valid (*conn));
  p4est_connectivity_complete (*conn);
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
octbingz2connectivity
(const char *filename, p4est_connectivity_t **conn)
{

  // load data from file
  octave_value tmp;
  octave_io_mode m = gz_read_mode;
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_load ("msh", tmp) == 0);

  Matrix p_matrix =
    tmp.scalar_map_value ().contents ("p").matrix_value ();

  Array<octave_idx_type> 
    t_matrix = tmp.scalar_map_value ().contents ("t").octave_idx_type_vector_value (false, true, true);
  //    t_matrix = tmp.scalar_map_value ().contents ("t").octave_idx_type_vector_value ();

  p4est_topidx_t num_vertices = p_matrix.cols (),
    num_trees = t_matrix.cols ();

  arrays2connectivity (p_matrix.fortran_vec (),
                       num_vertices,
                       t_matrix.fortran_vec (),
                       num_trees, conn);
};


void
tmesh::read_connectivity (const char *filename,
                          int source)
{
  if (rank == source)
    octbingz2connectivity (filename, &conn);

  conn = p4est_connectivity_bcast (conn, source, comm);
  p4est = p4est_new (comm, conn, 0, NULL, this);
};

void
tmesh::read_connectivity (const double *p,
                          const p4est_topidx_t num_vertices,
                          const p4est_topidx_t *t,
                          const p4est_topidx_t num_trees,
                          int source)
{
  if (rank == source)
    arrays2connectivity (p, num_vertices,
                         t, num_trees, &conn);

  conn = p4est_connectivity_bcast (conn, source, comm);
  p4est = p4est_new (comm, conn, 0, NULL, this);
};

void
tmesh::save (const char *filename)
{ p4est_save (filename, p4est, 0); };

void
tmesh::load (const char *filename)
{ p4est = p4est_load (filename, comm, 0, 0, this, &conn); };

void
tmesh::vtk_export (const char *filename)
{
  p4est_vtk_context_t *context =
    p4est_vtk_context_new (p4est, filename);
  assert (context != nullptr);
  p4est_vtk_context_set_scale (context, 1.0);
  p4est_vtk_context_set_continuous (context, 1);
  context = p4est_vtk_write_header (context);
  context =
    p4est_vtk_write_cell_dataf (context, 1, 1, 1, 0, 0, 0, context);
  assert (p4est_vtk_write_footer (context) == 0);
};

void
tmesh::octbin_export (const char * basename,
                      const std::vector<double> & f)
{
  assert (f.size () == num_global_nodes ());

  std::vector<double> p (2 * num_owned_nodes ());
  std::vector<double> f_loc (num_owned_nodes ());

  Array<octave_idx_type>
    oct_t (dim_vector (4, num_local_quadrants ()), 0);
  octave_idx_type *t = oct_t.fortran_vec ();

  std::array<tmesh::idx_t, 2> parents;
  std::array<int, 4> local_idx = {0, 1, 3, 2};

  octave_idx_type ij = 0;
  for (auto quadrant = begin_quadrant_sweep ();
       quadrant != end_quadrant_sweep ();
       ++quadrant)
    {
      ij = 0;
      for (auto ii : local_idx)
        {
          if ((! quadrant->is_hanging (ii))
              && (quadrant->t (ii) < num_owned_nodes ()))
            {
              p[2 * quadrant->t (ii) + 0] = quadrant->p (0, ii);
              p[2 * quadrant->t (ii) + 1] = quadrant->p (1, ii);
              f_loc[quadrant->t (ii)] = f[quadrant->gt (ii)];
              t[4 * quadrant->get_forest_quad_idx () + (ij++)] =
                quadrant->t (ii);
            }
          else if (! quadrant->is_hanging (ii))
            {
              p.push_back (quadrant->p (0, ii));
              p.push_back (quadrant->p (1, ii));
              f_loc.push_back (f[quadrant->gt (ii)]);
              t[4 * quadrant->get_forest_quad_idx () + (ij++)] =
                f_loc.size () - 1;
            }
          else if (quadrant->is_hanging (ii))
            {
              p.push_back (quadrant->p (0, ii));
              p.push_back (quadrant->p (1, ii));
              parents[0] = quadrant->gparent (0, ii);
              parents[1] = quadrant->gparent (1, ii);
              f_loc.push_back ((f[parents[0]] + f[parents[1]]) / 2.0);
              t[4 * quadrant->get_forest_quad_idx () + (ij++)] =
                f_loc.size () - 1;
            }
        }
    }


  Matrix oct_p (2, p.size () / 2, 0.0);
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
  sprintf (filename, "%s_%4.4d.octbin.gz", basename, rank);

  // Save to filename.
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_save ("msh", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);

};

void
tmesh::octbin_export_quadrant (const char * basename,
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
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_save ("msh", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);
};

void
tmesh::quadrant_iterator::reset ()
{
  if (data != nullptr)
    {
      p4est_t *p4 = data->the_tmesh->p4est;
      data->tree_idx         = p4->first_local_tree;
      data->tree_quad_idx    = 0;
      data->forest_quad_idx  = 0;

      data->is_ghost = false;
      data->qtq = -1;

      if (data->tree_idx != -1)
        {
          data->tree             =
            p4est_tree_array_index (p4->trees, data->tree_idx);
          data->tquadrants       = &(data->tree->quadrants);
          data->num_quadrants    =
            (p4est_locidx_t) data->tquadrants->elem_count;
        }
      else
        {
          data->tree             = nullptr;
          data->tquadrants       = nullptr;
          data->num_quadrants    = 0;
        }

      if (data->num_quadrants > 0)
        {
          auto tmp =
            p4est_quadrant_array_index (data->tquadrants,
                                        data->forest_quad_idx);
          data->update (data->tree_idx, tmp);
        }
      else
        data = nullptr;
    }
};

tmesh::quadrant_iterator
tmesh::begin_quadrant_sweep ()
{
  if (this->mesh == nullptr)
    this->update ();

  quadrant_iterator qi (&current_quadrant);
  qi.reset ();
  return qi;
};

void
tmesh::set_metrics_marker
(std::function<double (tmesh::quadrant_iterator)> estimator,
 double tol, int max_depth, int n_refine, int n_coarsen)
{
  this->metrics_max_depth = max_depth;

  int hxhat_hx = 0;

  for (auto quadrant = this->begin_quadrant_sweep ();
       quadrant != this->end_quadrant_sweep ();
       ++quadrant)
    {
      hxhat_hx = static_cast<int> (std::round (std::log2 (estimator (quadrant)
							  * std::sqrt (this->num_global_quadrants ()) / tol)));

      if (hxhat_hx >= 0)
	hxhat_hx = std::max (0, hxhat_hx - n_refine);
      else
	hxhat_hx = std::min (0, hxhat_hx + n_coarsen);
	
      quadrant->the_quadrant->p.user_int =
        std::min (std::max (-max_depth, hxhat_hx), max_depth);
    }

  return;
}

void
tmesh::refine (int recursive, int partforcoarsen, int balance)
{
  quadrant_iterator qi (&current_quadrant);
  qi.reset ();

  if (replace_fun == nullptr)
    p4est_refine (p4est, recursive, refine_callback, nullptr);
  else
    p4est_refine_ext (p4est, recursive, -1, refine_callback,
                      nullptr, replace_callback);

  if (balance)
    p4est_balance (p4est, P4EST_CONNECT_FACE, nullptr);

  p4est_partition (p4est, partforcoarsen, nullptr);

  if (! (lnodes == nullptr)) p4est_lnodes_destroy (lnodes);
  lnodes = nullptr;

  if (! (mesh == nullptr)) p4est_mesh_destroy (mesh);
  mesh = nullptr;

  if (! (ghost == nullptr)) p4est_ghost_destroy (ghost);
  ghost = nullptr;
}

void
tmesh::metrics_refine (idx_t max_elems)
{
  int recursive = 0;
  int partforcoarsen = 1;

  for (int i = 0; i < metrics_max_depth - 1; ++i)
    {
      coarsen (recursive, partforcoarsen, 0);
      refine (recursive, partforcoarsen, 0);

      // Prevent large meshes.
      if (max_elems > 0 && this->num_global_quadrants () >= max_elems)
        break;
    }

  coarsen (recursive, partforcoarsen, 0);
  refine (recursive, partforcoarsen, 1);
}

void
tmesh::coarsen (int recursive, int partforcoarsen, int balance)
{
  if (replace_fun == nullptr)
    p4est_coarsen (p4est, recursive, coarsen_callback, nullptr);
  else
    p4est_coarsen_ext (p4est, recursive, 0, coarsen_callback,
                       nullptr, replace_callback);

  if (balance)
    p4est_balance (p4est, P4EST_CONNECT_FACE, nullptr);

  p4est_partition (p4est, partforcoarsen, nullptr);

  if (! (lnodes == nullptr)) p4est_lnodes_destroy (lnodes);
  lnodes = nullptr;

  if (! (mesh == nullptr)) p4est_mesh_destroy (mesh);
  mesh = nullptr;

  if (! (ghost == nullptr)) p4est_ghost_destroy (ghost);
  ghost = nullptr;
};

void
tmesh::update ()
{
  ghost  = p4est_ghost_new  (p4est, P4EST_CONNECT_FULL);
  lnodes = p4est_lnodes_new (p4est, ghost, 1);
  mesh   = p4est_mesh_new   (p4est, ghost, P4EST_CONNECT_FACE);

  update_ghosts ();
}

void
tmesh::update_ghosts ()
{

  // deallocate previously built storage arrays
  if (! (this->mirror_data == nullptr))
    delete[] this->mirror_data;
  if (! (this->ghost_data  == nullptr))
    delete[] this->ghost_data;

  
  // Send mirror data.
  constexpr p4est_locidx_t chunk_len = 12;
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

      p4est_quadrant_t * q;

      // Loop over mirrors.
      for (p4est_locidx_t j = start; j < end; ++j)
        {
          q =
            p4est_quadrant_array_index (&ghost->mirrors,
                                        ghost->mirror_proc_mirrors[j]);

          p4est_tree_t * tree =
            p4est_tree_array_index (p4est->trees, q->p.which_tree);

          idx_t global_idx =
            p4est->global_first_quadrant[rank] +
            q->p.piggy3.local_num;

          quadrant_t current_mirror (this, q->p.which_tree, q);
          current_mirror.forest_quad_idx = q->p.piggy3.local_num;
          current_mirror.tree_quad_idx =
            current_mirror.forest_quad_idx - tree->quadrants_offset;
          current_mirror.update (q->p.which_tree, q);

          for (int node = 0; node < 4; ++node)
            mirror_data[mirror_end++] = current_mirror.gt (node);

          for (int node = 0; node < 4; ++node)
            if (current_mirror.is_hanging (node))
              {
                mirror_data[mirror_end++] =
                  current_mirror.gparent (0, node);
                mirror_data[mirror_end++] =
                  current_mirror.gparent (1, node);
              }
            else
              {
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
tmesh::user_int_replace (std::vector<int> old_user_int)
{
  std::vector<int> new_user_int;

  // Refinement.
  if (old_user_int.size () == 1)
    {
      new_user_int.resize (4);

      for (size_t i = 0; i < new_user_int.size (); ++i)
        new_user_int[i] = old_user_int[0] - 1;
    }
  // Coarsening.
  else if (old_user_int.size () == 4)
    {
      new_user_int.resize (1);

      new_user_int[0] = *std::max_element (old_user_int.begin (),
					   old_user_int.end ()) + 1;
    }

  return new_user_int;
}

int
tmesh::refine_callback (p4est_t* p4, p4est_topidx_t tt,
                        p4est_quadrant_t* qq)
{
  return (qq->p.user_int > 0);
};

int
tmesh::coarsen_callback (p4est_t* p4, p4est_topidx_t tt,
                         p4est_quadrant_t* qq [])
{
  return (qq[0]->p.user_int < 0
          && qq[1]->p.user_int < 0
          && qq[2]->p.user_int < 0
          && qq[3]->p.user_int < 0);
};

void
tmesh::replace_callback (p4est_t * p4,
                         p4est_topidx_t tt,
                         int num_outgoing,
                         p4est_quadrant_t * outgoing[],
                         int num_incoming,
                         p4est_quadrant_t * incoming[])
{
  tmesh *tm = reinterpret_cast<tmesh*> (p4->user_pointer);

  std::vector<int> old_user_int (num_outgoing);

  for (size_t i = 0; i < num_outgoing; ++i)
    old_user_int[i] = outgoing[i]->p.user_int;

  std::vector<int> new_user_int = tm->replace_fun (old_user_int);

  for (size_t i = 0; i < num_incoming; ++i)
    incoming[i]->p.user_int = new_user_int[i];

  return;
};
