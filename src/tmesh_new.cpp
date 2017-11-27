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

      data->num_quadrants = (p4est_locidx_t) data->tquadrants->elem_count;
    }

  auto tmp = p4est_quadrant_array_index (data->tquadrants, data->tree_quad_idx);
  data->update (data->tree_idx, tmp);
  
};

void 
tmesh::neighbor_iterator::operator++ ()
{
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  p4est_topidx_t which_tree;
  p4est_locidx_t which_quad;
  int nface, nrank;
  p4est_quadrant_t * neighbor = p4est_mesh_face_neighbor_next (data->face_neighbor,
                                                               &which_tree,
                                                               &which_quad,
                                                               &nface,
                                                               &nrank);
  
  if (neighbor != nullptr)
    {
      p4est_tree_t * tree = p4est_tree_array_index(data->the_tmesh->p4est->trees, which_tree);
      
      data->tree_idx = which_tree;
      
      if (data->face_neighbor->current_qtq < data->the_tmesh->num_local_quadrants())
        {
          data->forest_quad_idx = tree->quadrants_offset + which_quad;
          data->tree_quad_idx = which_quad;
        }
      else
        {
          data->forest_quad_idx = neighbor->p.piggy3.local_num +
                                  (data->the_tmesh->p4est->global_first_quadrant[nrank] -
                                   data->the_tmesh->p4est->global_first_quadrant[rank]);
          data->tree_quad_idx = neighbor->p.piggy3.local_num - tree->quadrants_offset;
        }
      data->forest_quad_idx = tree->quadrants_offset + which_quad;
      
      this->face_idx = nface;
      
      data->update(which_tree, neighbor);
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
      for (i = 0; i < 4; ++i)
        {
          tbuff[i] = ln->element_nodes[4 * forest_quad_idx + i];
          hbuff[i] = false;
          pbuff[i] = -1;
          pbuff[i+1] = -1;
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
  assert (pbuff[ip + in * 2] >= 0);

  return p4est_lnodes_global_index
    (the_tmesh->lnodes,
     static_cast<p4est_locidx_t>
     (tbuff[pbuff[ip + in * 2]]));
};

tmesh::idx_t
tmesh::quadrant_t::e (idx_t i)
{  
  assert (i < 4);
  idx_t retval = NOT_ON_BOUNDARY;
  p4est_quadrant_t node;
  p4est_quadrant_corner_node (this->the_quadrant, i, &node);
  
  if (node.y == 0)
    retval = 0;
  else if (node.y == P4EST_ROOT_LEN)
    retval = 1;
  else if (node.x == 0)
    retval = 2;
  else if (node.x == P4EST_ROOT_LEN)
    retval = 3;

  return retval;
};

tmesh::neighbor_iterator
tmesh::quadrant_t::begin_neighbor_sweep ()
{
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  
  if (! this->the_tmesh->mesh)
    {
      this->the_tmesh->update ();
    }
  
  p4est_mesh_face_neighbor_init (face_neighbor,
                                 this->the_tmesh->p4est,
                                 this->the_tmesh->ghost,
                                 this->the_tmesh->mesh,
                                 this->get_tree_idx(),
                                 this->the_quadrant);
  
  p4est_topidx_t which_tree;
  p4est_locidx_t which_quad;
  int nface, nrank;
  
  p4est_quadrant_t * neighbor = p4est_mesh_face_neighbor_next (face_neighbor,
                                                               &which_tree, &which_quad,
                                                               &nface, &nrank);
  
  current_neighbor = new quadrant_t(this->the_tmesh, which_tree, neighbor);

  p4est_tree_t * tree = p4est_tree_array_index(this->the_tmesh->p4est->trees, which_tree);
  
  current_neighbor->tree_idx = which_tree;
  
  if (face_neighbor->current_qtq < the_tmesh->num_local_quadrants())
    {
      current_neighbor->forest_quad_idx = tree->quadrants_offset + which_quad;
      current_neighbor->tree_quad_idx = which_quad;
    }
  else
    {
      current_neighbor->forest_quad_idx = neighbor->p.piggy3.local_num +
                                          (the_tmesh->p4est->global_first_quadrant[nrank] -
                                           the_tmesh->p4est->global_first_quadrant[rank]);
      current_neighbor->tree_quad_idx = neighbor->p.piggy3.local_num - tree->quadrants_offset;
    }
  current_neighbor->forest_quad_idx = tree->quadrants_offset + which_quad;
  current_neighbor->face_neighbor = this->face_neighbor;
  
  current_neighbor->update(which_tree, neighbor);
  
  neighbor_iterator ni (current_neighbor, nface);
  return ni;
}

tmesh::idx_t
tmesh::quadrant_t::t (tmesh::idx_t i)
{ return tbuff[i]; };

tmesh::idx_t
tmesh::quadrant_t::gt (tmesh::idx_t i)
{
  return p4est_lnodes_global_index
    (the_tmesh->lnodes,
     static_cast<p4est_locidx_t> (tbuff[i]));
};

bool
tmesh::quadrant_t::is_hanging (tmesh::idx_t i)
{ return hbuff[i]; };


tmesh::~tmesh ()
{
  p4est_destroy (this->p4est);
  p4est_connectivity_destroy (this->conn);
  if (!(this->lnodes == nullptr)) p4est_lnodes_destroy (this->lnodes);
  if (!(this->mesh   == nullptr)) p4est_mesh_destroy   (this->mesh);
  if (!(this->ghost  == nullptr)) p4est_ghost_destroy  (this->ghost);
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
tmesh::octbin_export (const char * basename,
                      const std::vector<double> & f,
                      MPI_Comm comm)
{
  assert (f.size () == num_global_nodes ());
    
  std::vector<double> p (2 * num_owned_nodes ());
  std::vector<double> f_loc (num_owned_nodes ());  

  Array<octave_idx_type> oct_t (dim_vector (4, num_local_quadrants ()), 0);
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
          if ((! quadrant->is_hanging (ii)) && (quadrant->t (ii) < num_owned_nodes ()))
            {            
              p[2 * quadrant->t (ii) + 0] = quadrant->p (0, ii);
              p[2 * quadrant->t (ii) + 1] = quadrant->p (1, ii);            
              f_loc[quadrant->t (ii)] = f[quadrant->gt (ii)];
              t[4 * quadrant->get_forest_quad_idx () + (ij++)] = quadrant->t (ii);
            }
          else if (! quadrant->is_hanging (ii))
            {
              p.push_back (quadrant->p (0, ii));
              p.push_back (quadrant->p (1, ii));
              f_loc.push_back (f[quadrant->gt (ii)]);
              t[4 * quadrant->get_forest_quad_idx () + (ij++)] = f_loc.size () - 1;
            }
          else if (quadrant->is_hanging (ii))
            {
              p.push_back (quadrant->p (0, ii));
              p.push_back (quadrant->p (1, ii));
              parents[0] = quadrant->gparent (0, ii);
              parents[1] = quadrant->gparent (1, ii);
              f_loc.push_back ((f[parents[0]] + f[parents[1]]) / 2.0);
              t[4 * quadrant->get_forest_quad_idx () + (ij++)] = f_loc.size () - 1;
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
  int rank;
  MPI_Comm_rank (comm, &rank);

  char filename[255] = "";
  sprintf(filename, "%s_%4.4d.octbin.gz", basename, rank);
  
  // Save to filename.
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_save ("msh", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);

};

void
tmesh::quadrant_iterator::reset ()
{
  p4est_t *p4 = data->the_tmesh->p4est;
  data->tree_idx         = p4->first_local_tree;
  data->tree_quad_idx    = 0;
  data->forest_quad_idx  = 0;

  data->tree             = p4est_tree_array_index (p4->trees, data->tree_idx);
  data->tquadrants       = &(data->tree->quadrants);
  data->num_quadrants    = (p4est_locidx_t) data->tquadrants->elem_count;

  if (data->num_quadrants > 0)
    {
      auto tmp = p4est_quadrant_array_index (data->tquadrants, data->forest_quad_idx);
      data->update (data->tree_idx, tmp);
    }
  else
    data = nullptr;  
};

tmesh::quadrant_iterator
tmesh::begin_quadrant_sweep ()
{
  quadrant_iterator qi (&current_quadrant);
  qi.reset ();
  return qi;
};

void
tmesh::refine (int recursive, int partforcoarsen)
{
  p4est_refine (p4est, recursive, refine_callback, nullptr);
  p4est_balance (p4est, P4EST_CONNECT_FULL, nullptr);
  p4est_partition (p4est, partforcoarsen, nullptr);
  
  if (! (lnodes == nullptr)) p4est_lnodes_destroy (lnodes);
  lnodes = nullptr;

  if (! (mesh == nullptr)) p4est_mesh_destroy (mesh);
  mesh = nullptr;
}

void
tmesh::coarsen (int recursive, int partforcoarsen)
{
  p4est_coarsen (p4est, recursive, coarsen_callback, nullptr);
  p4est_balance (p4est, P4EST_CONNECT_FULL, nullptr);
  p4est_partition (p4est, partforcoarsen, nullptr);

  if (! (lnodes == nullptr)) p4est_lnodes_destroy (lnodes);
  lnodes = nullptr;  

  if (! (mesh == nullptr)) p4est_mesh_destroy (mesh);
  mesh = nullptr;
};

void
tmesh::update ()
{
  ghost  = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);
  lnodes = p4est_lnodes_new (p4est, ghost, 1);
  mesh   = p4est_mesh_new (p4est, ghost, P4EST_CONNECT_FULL);
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
