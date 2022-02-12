#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <map>
#include <vector>
#include <nanoshaper.h>

using int_coord_t = unsigned long long int;

struct
crossings_t {

   /*struct
   vector3d {
     double x;
     double y;
     double z;
   };*/

   double start[3],  end; //start and end point of the choosen direction
   //bool   startv, endv; //??

   //std::vector<double>   params; 
   std::vector<bool>     flags; //in or out (in=1, out=0) ex: if flags[1]==1 -> [inters[1].first ; inters[2].first] is IN  
   //std::vector<vector3d> normals; //vector of normals 
   
   std::vector<std::pair<double,double*> > inters; //intersections and normals
   
   crossings_t (double start_x, double start_y, double start_z, double end_, NS::NanoShaper & ns_)
    : end(end_)
    {
      start[0] = start_x;
      start[1] = start_y;
      start[2] = start_z;
      computeIntersections (ns_);
    };
   
   void
   computeIntersections (NS::NanoShaper & ns); //it fills inters vector
   
   double
   is_inside_molecule(double point); //given a point coord it returns 1 if the point is inside the molecule, 0 otherwise

};

/*
struct
ray_cache_t 
   {
   std::map<int_coord_t, std::map<int_coord_t, crossings_t>> rays;
   
   crossings_t &
   operator() (int_coord_t x0, int_coord_t x1) 
   {
     static crossings_t dummy;
     crossings_t & tmp = dummy;
     auto it0 = rays.find (x0);
     
     if (it0 != rays.end ()) //if x0 is present in rays
     {
       auto it1 = it0->second.find (x1); //it1 is the map with x1 as first elem
       if (it1 != it0->second.end ())  //if it1 exists
         tmp = it1->second; //temp is crossings_t of it1
     }
      
     else //x0 not in rays 
     {
       // create a new ray
       crossings_t tmp;
       rays[x0][x1] = tmp; //create a map element with first elem x0, then the second element is a map with x1 as first, tmp as second elem
       tmp = rays[x0][x1];
     }
     return tmp;
   }
   
};*/

#endif //RAYTRACER_H
