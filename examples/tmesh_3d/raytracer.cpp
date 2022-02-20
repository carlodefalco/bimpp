#include "raytracer.h"

void
crossings_t::compute_intersections ()
{
  unsigned y_direction = 1;
  
  double start_ray[3] = {point[0], start, point[1]};
  bool compute_normals = false;
  
  if (point[0] < start || point[1] < start || point[0] > end || point[1] > end) //if I'm certainly out of the molecule
    return;
    
  ns.castAxisOrientedRay (start_ray, end, inters, y_direction, compute_normals);

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
crossings_t::is_inside_molecule (double y)
{
  unsigned i = 0;
  if (flags.size () == 0 || y < inters[i].first)
     return 0; //if there are no inters or y_the coord is before the first intersection, the point is outside.
    
  while (i < flags.size () && y > inters[i].first) //go on until the inters is passed
     i++;
     
  return flags[i-1];

}
