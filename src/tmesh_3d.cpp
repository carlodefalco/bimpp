/*! \file tmesh_3d.cpp
  \brief Interface for p8est library
*/

#include <tmesh_3d.h>
#include <array>

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

