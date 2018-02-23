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
static const int    zero = 0;      /**< Constant zero. */
static const int    ones = P8EST_CHILDREN - 1;  /**< One bit per dimension. */
static const int   *corner_to_hanging[P8EST_CHILDREN];
static const int    corner_num_hanging[P8EST_CHILDREN] = { 1, 2, 2, 4, 2, 4, 4, 1 };

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

tmesh_3d::~tmesh_3d ()
{
    p8est_destroy (this->p8est);
    p8est_connectivity_destroy (this->conn);
    if (! (this->lnodes == nullptr)) p8est_lnodes_destroy (this->lnodes);
    if (! (this->mesh   == nullptr)) p8est_mesh_destroy   (this->mesh);
    if (! (this->ghost  == nullptr)) p8est_ghost_destroy  (this->ghost);
    
    if (! (this->mirror_data == nullptr)) delete[] this->mirror_data;
    if (! (this->ghost_data  == nullptr)) delete[] this->ghost_data;
};

