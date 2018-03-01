/*! \file tmesh_3d.cpp
  \brief Interface for p8est library
*/

#include <tmesh_3d.h>
#include <array>

double
tmesh_3d::octant_t::p (tmesh_3d::idx_t ii, tmesh_3d::idx_t jj) 
{
  double retval = vxyz[P8EST_DIM*jj+ii];
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
static const int    zero = 0;      /**< Constant zero. */
static const int    ones = P8EST_CHILDREN - 1;  /**< One bit per dimension. */
static const int   *corner_to_hanging[P8EST_CHILDREN];
static const int    corner_num_hanging[P8EST_CHILDREN] = { 1, 2, 2, 4, 2, 4, 4, 1 };

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
{//CESARE(this is just a name translation, check for what neighbor means)
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
  int hanging_corner[P8EST_CHILDREN];
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
  
  for (i = 0; i < P8EST_CHILDREN; ++i)
    {
      p8est_quadrant_corner_node (this->the_octant, i, &node);
      p8est_qcoord_to_vertex (this->the_tmesh->conn, tree_idx,
                              node.x, node.y, node.z, &(vxyz[P8EST_DIM * i]));
    }

  if (ln != nullptr)
    {
      // Non-ghost elements.
      if (! is_ghost)
        {
          for (i = 0; i < P8EST_CHILDREN; ++i)
            {
              tbuff[i] = ln->element_nodes[P8EST_CHILDREN * forest_oct_idx + i];
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
            for (i = 0; i < P8EST_CHILDREN; ++i)
              if (hanging_corner[i] >= 0)
                {
                  hbuff[i] = true;
                  c = hanging_corner[i];
                  num_parents = corner_num_hanging[i ^ c];
                  base_corner = corner_to_hanging[i ^ c];
                  for (j = 0; j < num_parents; ++j)
                    pbuff[j + 4 * i] = base_corner[j] ^ c;//CESARE(???)
                }
        }
      // Ghost elements.
      else
        {
          p4est_locidx_t idx = this->qtq -
            the_tmesh->num_local_octants ();
          
          for (i = 0; i < P8EST_CHILDREN; ++i)
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

std::vector<int>
tmesh_3d::userint_replace (std::vector<int> old_userint)
{
  std::vector<int> new_userint;
  
  // Refinement.
  if (old_userint.size () == 1)
    {
      new_userint.resize (P8EST_CHILDREN);
      
      for (size_t i = 0; i < new_userint.size (); ++i)
        new_userint[i] = old_userint[0] - 1;
    }
  // Coarsening.
  else if (old_userint.size () == P8EST_CHILDREN)
    {
      new_userint.resize (1);
      
      new_userint[0] = *std::min_element (old_userint.begin (),
                                          old_userint.end ()) + 1;
    }
  
  return new_userint;
}

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

