#include "pb_class.h"
#include "GetPot"


#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>


#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>


#include <p8est.h>

void
poisson_boltzmann::read_csv ()
{
  std::string filename = csvfilename;
  std::ifstream in;
  std::string buf;

  in.open (filename);
  if (! in) std::cerr << "could not open file " << filename << std::endl;
  try
    {
      while (std::getline (in, buf, ','))
        {
          double x = stod (buf);
          ions.push_back (ion ({x, .0, .0}, .0));
        
          std::getline (in, buf, ',');
          x = stod (buf);
          ions.back ().center[1] = x;

          std::getline (in, buf, ',');
          x = stod (buf);
          ions.back ().center[2] = x;

          std::getline (in, buf);
          x = stod (buf);
          ions.back ().radius = x;
        }
    }
  catch(std::exception &e)
    {
      std::cerr << "encountered exception " << e.what () << " while reading file " << filename << std::endl;
      exit (1);
    }

  auto comp = [] (const ion &i1, const ion &i2) -> bool { return i1.radius < i2.radius; };
  double maxradius = std::max_element (ions.begin (), ions.end (), comp)->radius;

  ll = 0; rr = 0; 
  auto it = [this] (const ion &i1)
    {
      for (int kk = 0; kk < 3; ++kk)
        {
          if (i1.center[kk] > this->rr)
            this->rr = i1.center[kk];
          else if (i1.center[kk] < this->ll)
            this->ll = i1.center[kk];
        }
    };
  std::for_each (ions.begin (), ions.end (), it);
  ll -= 3*maxradius;
  rr += 3*maxradius;
  
  simple_conn_p = {ll, ll, ll, rr, ll, ll, ll, rr, ll, rr, rr, ll,
                   ll, ll, rr, rr, ll, rr, ll, rr, rr, rr, rr, rr};

  simple_conn_t = {1, 2, 3, 4, 5, 6, 7, 8, 1};
  tmsh.read_connectivity (simple_conn_p.data (), simple_conn_num_vertices,
                          simple_conn_t.data (), simple_conn_num_trees);

}


double
poisson_boltzmann::levelsetfun (double x, double y, double z)
{
  double dist = 0.0; 
  for (const ion& i : ions)
    {
      
      dist += std::exp (decay * ((std::pow (x - i.center[0], 2) +
                                  std::pow (y - i.center[1], 2) +
                                  std::pow (z - i.center[2], 2)) /
                                 std::pow (i.radius, 2) - 1.0));

      if (dist > 1.5)
        break;
      
    }
  return dist;
}


void
poisson_boltzmann::parse_options (int argc, char **argv)
{
  GetPot g (argc, argv);
  csvfilename = g ("csvfilename", "benzene.csv");
  p4estfilename = g ("p4estfilename", "poisson_boltzmann_p4est");
  markerfilename = g ("markerfilename", "poisson_boltzmann_marker");
  minlevel = g ("minlevel", minlevel);
  maxlevel = g ("maxlevel", maxlevel);
  decay    = g ("decay", decay);
}

void
poisson_boltzmann::init_tmesh ()
{
  for (auto i = 0; i < minlevel; ++i)
    {
      tmsh.set_refine_marker (uniform_refinement);       
      tmsh.refine (0, 1);      
    }
}


void
poisson_boltzmann::refine_surface ()
{
  for (int kk = 0; kk < (maxlevel - minlevel); ++kk)
    {
      // REFINEMENT
      {

        distributed_vector rcoeff (tmsh.num_owned_nodes ());
  
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant)
          {

            for (int ii = 0; ii < 8; ++ii)
              {
                if (! quadrant->is_hanging (ii))
                  rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                                           quadrant->p (1, ii),
                                                           quadrant->p (2, ii));
                else
                  for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                    rcoeff[quadrant->gparent (jj, ii)] += 0.;
              }
          }
        //bim3a_solution_with_ghosts (tmsh, rcoeff, replace_op);
  
        auto refinement = [&rcoeff,this]
          (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 1.0;
            double min = 100.0 * this->ions.size ();
            double max = 0.0;
            double tmp = 0.0;

            if (currentlevel >= this->maxlevel)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {

                    if (! q->is_hanging (ii))
                      {
                        tmp = rcoeff[q->gt (ii)];

                        if (tmp > max) max = tmp;
                        if (tmp < min) min = tmp;
                      }
     
                  }
                if (max > 1.0 && min < 1.0)
                  retval = this->maxlevel - currentlevel;
              }

            return (retval);
          };

        tmsh.set_refine_marker (refinement);
        tmsh.refine (0, 1);      
      }

      // COARSENING
      {
        distributed_vector rcoeff (tmsh.num_owned_nodes ());
  
        for (auto quadrant = tmsh.begin_quadrant_sweep ();
             quadrant != tmsh.end_quadrant_sweep ();
             ++quadrant)
          {

            for (int ii = 0; ii < 8; ++ii)
              {
                if (! quadrant->is_hanging (ii))
                  rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p(0, ii),
                                                           quadrant->p(1, ii),
                                                           quadrant->p(2, ii));
                else
                  for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                    rcoeff[quadrant->gparent (jj, ii)] += 0;
              }
          }
        //bim3a_solution_with_ghosts (tmsh, rcoeff, replace_op);

        auto coarsening = [&rcoeff,this]
          (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 0;
            double min = 100.0 * this->ions.size ();
            double max = 0.0;
            double tmp = 0.0;

            if (currentlevel <= this->minlevel)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {

                    if (! q->is_hanging (ii))
                      {
                        tmp = rcoeff[q->gt (ii)];
                        
                        if (tmp > max) max = tmp;
                        if (tmp < min) min = tmp;
                      }
     
                  }

                if (min > 1.0 || max < 1.0)
                  retval = currentlevel - this->minlevel;
              }

            return (retval);
          };
        
        tmsh.set_coarsen_marker (coarsening);
        tmsh.coarsen (0, 1);      
      }
    }
}

void
poisson_boltzmann::export_marked_tmesh ()
{
    distributed_vector rcoeff (tmsh.num_owned_nodes ());
  
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
      {
            
        for (int ii = 0; ii < 8; ++ii)
          {
            if (! quadrant->is_hanging (ii))
              rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p(0, ii),
                                                       quadrant->p(1, ii),
                                                       quadrant->p(2, ii));
            else
              for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                rcoeff[quadrant->gparent (jj, ii)] += 0.;
          }
      }
    bim3a_solution_with_ghosts (tmsh, rcoeff, replace_op);
    tmsh.octbin_export (markerfilename.c_str (), rcoeff);
}

void
poisson_boltzmann::export_p4est ()
{
  tmsh.save (p4estfilename.c_str ());
}
  
