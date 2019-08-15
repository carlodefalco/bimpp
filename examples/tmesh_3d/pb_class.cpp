#include "pb_class.h"
#include "GetPot"


#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>
#include <mumps_class.h>


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

          std::getline (in, buf, ',');
          x = stod (buf);
          ions.back ().radius = x;

          std::getline (in, buf);
          x = stod (buf);
          ions.back ().charge = x;
        }
    }
  catch(std::exception &e)
    {
      std::cerr << "encountered exception "
                << e.what ()
                << " while reading file "
                << filename << std::endl;
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
  markerfilename = g ("markerfilename", "poisson_boltzmann_marker_0");
  lsfilename = g ("lsfilename", "poisson_boltzmann_levelset_0");
  minlevel = g ("minlevel", minlevel);
  maxlevel = g ("maxlevel", maxlevel);
  decay    = g ("decay", decay);
  e_in     = g ("e_in", e_in);
  e_out    = g ("e_out", e_out);
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

bool
poisson_boltzmann::is_in (const poisson_boltzmann::ion& i,
                          tmesh_3d::quadrant_iterator q)
{
  double tol =  p4esttol * (rr-ll);
  bool retval = false;
  double l, r, t, b, f, bk;

  l = q->p (0, 0);
  r = q->p (0, 0);

  f  = q->p (1, 0);
  bk = q->p (1, 0);

  b = q->p (2, 0);
  t = q->p (2, 0);


  for (int ii = 1; ii < 8; ++ii)
    {
      l  = q->p (0, ii)  < l  ? q->p (0, ii) : l;
      r  = q->p (0, ii)  > r  ? q->p (0, ii) : r;
      f  = q->p (1, ii)  < f  ? q->p (1, ii) : f;
      bk = q->p (1, ii)  > bk ? q->p (1, ii) : bk;
      b  = q->p (2, ii)  < b  ? q->p (2, ii) : b;
      t  = q->p (2, ii)  > t  ? q->p (2, ii) : t;

    }

  retval =           (i.center[0] > l - tol) && (i.center[0] <= r  - tol);
  retval = retval && (i.center[1] > f - tol) && (i.center[1] <= bk - tol);
  retval = retval && (i.center[2] > b - tol) && (i.center[2] <= t  - tol);

  return retval;
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
                  rcoeff[quadrant->gt (ii)] =
                    levelsetfun (quadrant->p (0, ii),
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
                else
                  for (const ion& i : ions)
                    if (is_in (i, q))
                      {
                        retval = this->maxlevel - currentlevel;
                        break;
                      }
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

                for (const ion& i : ions)
                  if (is_in (i, q))
                    {
                      retval = 0;
                      break;
                    }
              }

            return (retval);
          };

        tmsh.set_coarsen_marker (coarsening);
        tmsh.coarsen (0, 1);
      }
    }
}

void
poisson_boltzmann::create_markers ()
{
  this->marker.assign (this->tmsh.num_local_quadrants (), 0.0);
  this->rho_fixed.assign (this->tmsh.num_local_quadrants (), 1.0);
  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
       quadrant != this->tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (const ion& i : ions)
        if (is_in (i, quadrant))
          {
            this->rho_fixed[quadrant->get_forest_quad_idx ()] = 1.0e11;
            break;
          }

      int num_int_nodes = 0;
      int num_hanging = 0;
      for (int ii = 0; ii < 8; ++ii)
        {
          if (! quadrant->is_hanging (ii))
            {
              if (this->levelsetfun (quadrant->p (0, ii),
                                     quadrant->p (1, ii),
                                     quadrant->p (2, ii)) > 1.0)
                ++num_int_nodes;
            }
          else
            ++num_hanging;
        }
      if (num_int_nodes == 0)
        this->marker[quadrant->get_forest_quad_idx ()] = 1.0;
      else if (num_int_nodes < (8 - num_hanging))
        this->marker[quadrant->get_forest_quad_idx ()] = 1.0/2.0;
    }
}

void
poisson_boltzmann::export_ls_tmesh ()
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
    tmsh.octbin_export (lsfilename.c_str (), rcoeff);
}

void
poisson_boltzmann::export_marked_tmesh ()
{
  tmsh.octbin_export_quadrant (markerfilename.c_str (), marker);
}

void
poisson_boltzmann::export_p4est ()
{
  tmsh.save (p4estfilename.c_str ());
}

void
poisson_boltzmann::compute_electric_potential ()
{
  epsilon.assign (tmsh.num_local_quadrants (), e_in);
  for (auto epsp = epsilon.begin (), mp = marker.begin ();
       epsp != epsilon.end () || mp != marker.end ();
       ++epsp, ++mp)
    if ((*mp) == 0.0)
      (*epsp) = e_in;
    else
      (*epsp) = e_out;

  reaction.assign (tmsh.num_local_quadrants (), 0.0);
  for (auto rp = reaction.begin (), mp = marker.begin ();
       rp != reaction.end () || mp != marker.end ();
       ++rp, ++mp)
    if ((*mp) != 0.0)
      (*rp) = 1.22e-06;

  tmsh.octbin_export_quadrant ("epsilon_0", epsilon);
  tmsh.octbin_export_quadrant ("rho_0", rho_fixed);
  tmsh.octbin_export_quadrant ("reaction_0", reaction);

  sparse_matrix A;
  A.resize (tmsh.num_global_nodes ());
  distributed_vector  rhs (tmsh.num_global_nodes ());

  distributed_vector  psi (tmsh.num_global_nodes ());
  psi.get_owned_data ().assign (psi.get_owned_data ().size (), 0.0);
  bim3a_solution_with_ghosts (tmsh, psi);
  distributed_vector ones (tmsh.num_global_nodes ());
  ones.get_owned_data ().assign (ones.get_owned_data ().size (), 1.0);
  bim3a_solution_with_ghosts (tmsh, ones, replace_op);
  bim3a_advection_diffusion (tmsh, epsilon, psi, A);
  bim3a_reaction (tmsh, reaction, ones, A);


  bim3a_rhs (tmsh, rho_fixed, ones, rhs);
  tmsh.octbin_export ("rhs_0", rhs);
  
  mumps mumps_solver;

  std::vector<double> vals;
  std::vector<int> irow, jcol;

  A.aij (vals, irow, jcol, mumps_solver.get_index_base ());

  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);
  mumps_solver.set_rhs_distributed (rhs);

  std::cout << "mumps_solver.analyze () = "
            << mumps_solver.analyze ()
            << std::endl;
  std::cout << "mumps_solver.factorize () = "
            << mumps_solver.factorize ()
            << std::endl;
  std::cout << "mumps_solver.solve () = "
            << mumps_solver.solve ()
            << std::endl;

  distributed_vector phi = mumps_solver.get_distributed_solution ();
  bim3a_solution_with_ghosts (tmsh, phi);

  tmsh.octbin_export ("phi_0", phi);

  mumps_solver.cleanup ();
}
