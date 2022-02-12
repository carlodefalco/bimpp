#include "raytracer.h"

void
crossings_t::computeIntersections (NS::NanoShaper & ns)
{

  //set ray direction:
  unsigned y_direction = 1; 
  ns.setDirection(y_direction);
  
  bool computeNormals = false;
  ns.castAxisOrientedRay (start, end, inters, y_direction, computeNormals);
  
  if (inters.size() != 0)
  {
     for (unsigned i = 0; i < ((inters.size()) - 1); i++)
     {
        flags.push_back (1);
        flags.push_back (0);
        i++;
     }
  }

} 

double
crossings_t::is_inside_molecule(double y_coord)
{

  unsigned i = 0;
  if (flags.size () == 0 || y_coord < inters[i].first)
     return 0; //if there are no inters or y_the coord is before the first intersection, the point is outside.
    
  while (i < flags.size () && y_coord > inters[i].first) //go on until the inters is passed
     i++;
     
  return flags[i-1];

}
