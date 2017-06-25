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
tmesh::refine (int recursive, int partforcoarsen)
{
  p4est_refine (p4est, recursive, refine_callback, nullptr);
  p4est_balance (p4est, P4EST_CONNECT_FACE, nullptr);
  p4est_partition (p4est, partforcoarsen, nullptr);
}


tmesh::~tmesh ()
{
  p4est_destroy (this->p4est);
  p4est_connectivity_destroy (this->conn);
};

double
tmesh::quadrant_t::p (tmesh::idx_t ii, tmesh::idx_t jj) 
{
  double retval = vxyz[3*jj+ii];
  return (retval);
};



int
tmesh::refine_callback (p4est_t* p4, p4est_topidx_t tt, p4est_quadrant_t* qq)
{
  tmesh *tm = reinterpret_cast<tmesh*> (p4->user_pointer);
  tm->update_quadrant (tt, qq);
  quadrant_iterator qi (&(tm->current_quadrant));
  return tm->refine_marker (qi);
};

int
tmesh::coarsen_callback (p4est_t* p4, p4est_topidx_t tt, p4est_quadrant_t* qq[])
{
  tmesh *tm = reinterpret_cast<tmesh*> (p4->user_pointer);
};

void
tmesh::quadrant_t::update (p4est_topidx_t tree,
                           p4est_quadrant_t *q)
{
  p4est_quadrant_t node;
  idx_t i;
  this->the_tree = tree;
  this->the_quadrant = q;
  for (i = 0; i < 4; ++i)
    {
      p4est_quadrant_corner_node (this->the_quadrant, i, &node);
      p4est_qcoord_to_vertex (this->the_tmesh->conn, the_tree, node.x, node.y, &(vxyz[3 * i]));
    }
};

void
tmesh::update_quadrant (p4est_topidx_t tree,
                        p4est_quadrant_t *q)
{ current_quadrant.update (tree, q); };

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
