#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <map>
#include <vector>
#include <array>
#include <nanoshaper.h>

#include "json.hpp" 

using int_coord_t = unsigned long long int;

struct
crossings_t {

   static double start, end; //start and end point of the choosen direction (l_c[1], r_c[1])
   static NS::NanoShaper ns;
   
   double point[2]; //point x and z coords: the ones that prescribe the ray
 
   std::vector<bool> flags; //in or out (in=1, out=0) ex: if flags[1]==1 -> [inters[1].first ; inters[2].first] is IN  
   std::vector<std::pair<double,double*>> inters; //intersections and normals
   
   void
   compute_intersections (); //it fills inters vector
   
   double
   is_inside_molecule (double y);

};

struct
map_compare
{
  bool operator ()(const std::array<double, 2> & arr1, const std::array<double, 2> & arr2) const
  {  
    if (arr1.at(0) < arr2.at(0))
      return 1;
      
    else if (arr1.at(0) == arr2.at(0))
    {
      if (arr1.at(1) < arr2.at(1))
        return 1;
      else 
        return 0;
    }
    
    else 
      return 0;
  }
};

struct
ray_cache_t 
{
   std::map<std::array<double, 2> , crossings_t, map_compare> rays; //x-coord and z-coord, correspondent ray
   
   static int_coord_t count_cache;
   static int_coord_t count_new;
   
   crossings_t 
   operator() (double x0, double x1) 
   {
     static crossings_t dummy;
     crossings_t & cr_t = dummy;
     
     static std::array<double, 2> start_point;
     start_point = {x0, x1};
     auto it0 = rays.find (start_point);
     
     if (it0 != rays.end ())
     {
       cr_t = it0->second;
       count_cache++;
     }
      
     else
     {
       // create a new ray
       crossings_t tmp;
       tmp.point[0] = x0;
       tmp.point[1] = x1;
       tmp.compute_intersections ();
       rays[start_point] = tmp;         
       cr_t = rays[start_point];
       count_new++;
     }

     return cr_t;
   }

};
#endif //RAYTRACER_H
