#include "pb_class.h"
#include "GetPot"

#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>
#include <mumps_class.h>
#include <lis_class.h>

#include <cmath>
#include <cstdio>
#include <fstream>

#include <p8est.h>

void
poisson_boltzmann::create_cubic_mesh ()
{
  
  auto comp = [] (const NS::Atom &a1, const NS::Atom &a2) -> bool { return a1.radius < a2.radius; }; 
  double maxradius = std::max_element (atoms.begin (), atoms.end (), comp)->radius; 
  
  ll = 0; rr = 0;  
  auto it = [this] (const NS::Atom &a1)
    {
      for (int kk = 0; kk < 3; ++kk)
        {
          if (a1.pos[kk] > this->rr)
            this->rr = a1.pos[kk]; //trovo il punto più a dx
          else if (a1.pos[kk] < this->ll)
            this->ll = a1.pos[kk]; //trovo il punto più a sx
        }
    };
  std::for_each (atoms.begin (), atoms.end (), it);
  ll -= 3*maxradius;
  rr += 3*maxradius;

  simple_conn_p = {ll, ll, ll, rr, ll, ll, ll, rr, ll, rr, rr, ll,
                   ll, ll, rr, rr, ll, rr, ll, rr, rr, rr, rr, rr}; 

  simple_conn_t = {1, 2, 3, 4, 5, 6, 7, 8, 1};
  tmsh.read_connectivity (simple_conn_p.data (), simple_conn_num_vertices,
                          simple_conn_t.data (), simple_conn_num_trees);

}


void
poisson_boltzmann::create_mesh ()
{
  
  auto comp = [] (const NS::Atom &a1, const NS::Atom &a2) -> bool { return a1.radius < a2.radius; }; 
  double maxradius = std::max_element (atoms.begin (), atoms.end (), comp)->radius; 
  
  l_c[0] = 0.0; l_c[1] = 0.0; l_c[2] = 0.0; 
  r_c[0] = 0.0; r_c[1] = 0.0; r_c[1] = 0.0;
  auto it = [this] (const NS::Atom &a1)
    {
      for (int kk = 0; kk < 3; ++kk)
        {
          if (a1.pos[kk] > this->r_c[kk])
            this->r_c[kk] = a1.pos[kk]; 
          else if (a1.pos[kk] < this->l_c[kk])
            this->l_c[kk] = a1.pos[kk]; //trovo il punto più a sx
        }
    };
  std::for_each (atoms.begin (), atoms.end (), it);
  l_c[0] -= 3*maxradius; l_c[1] -= 3*maxradius; l_c[2] -= 3*maxradius;
  r_c[0] += 3*maxradius; r_c[1] += 3*maxradius; r_c[2] += 3*maxradius;

  simple_conn_p = {l_c[0], l_c[1], l_c[2], r_c[0], l_c[1], l_c[2], l_c[0], r_c[1], l_c[2], r_c[0], r_c[1], l_c[2],
                   l_c[0], l_c[1], r_c[2], r_c[0], l_c[1], r_c[2], l_c[0], r_c[1], r_c[2], r_c[0], r_c[1], r_c[2]}; 

  simple_conn_t = {1, 2, 3, 4, 5, 6, 7, 8, 1};
  tmsh.read_connectivity (simple_conn_p.data (), simple_conn_num_vertices,
                          simple_conn_t.data (), simple_conn_num_trees);

}

double
poisson_boltzmann::levelsetfun (double x, double y, double z)
{
  double dist = 0.0;
  for (const NS::Atom& i : atoms)
    {

      dist += std::exp (decay * ((std::pow (x - i.pos[0], 2) +
                                  std::pow (y - i.pos[1], 2) +
                                  std::pow (z - i.pos[2], 2)) /
                                  (i.radius*i.radius) - 1.0));
                                 
      if (dist > 1.5)
        break;

    }

  return dist;
}


int
poisson_boltzmann::parse_options (int argc, char **argv)
{
  GetPot g (argc, argv);
  optionsfile = g ("filepot", "../options.pot");
  pqrfile = g ("pqrfilename", "1CCM.pqr");
  
  // controlla che i file esistano  
  std::ifstream check1(optionsfile);
  if(!check1)
  {
  	std::cerr << "Cannot find the options file" << std::endl;
  	return 1;
  }
  
  std::ifstream check2(pqrfile);
  if(!check2)
  {
  	std::cerr << "Cannot find the pqr file" << std::endl;
  	return 1;
  }
  
  // leggi dal file le opzioni
  GetPot g2 (optionsfile.c_str ());
  
  const std::string mesh_options = "mesh/";
  maxlevel = g2 ((mesh_options + "maxlevel").c_str (),  6);
  minlevel = g2 ((mesh_options + "minlevel").c_str (),  4);
  
  const std::string model_options = "model/";
  linearized = g2 ((model_options + "linear_solver").c_str (),  1);
  ionic_strength = g2 ((model_options + "ionic_strength").c_str (),  0.145);
  e_in = g2 ((model_options + "molecular_dielectric_constant").c_str (),  2.);
  e_out = g2 ((model_options + "solvent_dielectric_constant").c_str (),  78.54);
  decay    = g2 ((model_options + "decay").c_str (), -1.5);
  
  const std::string alg_options = "algorithm/";
  linear_solver_name = g2 ((alg_options + "linear_solver").c_str (),  "mumps");
  linear_solver_options = g2 ((alg_options + "solver_options").c_str (),  "");

  const std::string out_options = "output/";
  p4estfilename = g2 ((out_options + "p4estfilename").c_str (), "poisson_boltzmann_p4est");
  markerfilename = g2 ((out_options + "markerfilename").c_str (), "poisson_boltzmann_marker_0");
  lsfilename = g2 ((out_options + "lsfilename").c_str (), "poisson_boltzmann_levelset_0");
  
  return 0;
}

void 
poisson_boltzmann::print_options ()
{
  std::cout << "\nChoosen options: " << std::endl;
  std::cout << "maxlevel = " << maxlevel <<  "; minlevel = " << minlevel << std::endl;
  std::cout << "Linearized model = " << linearized << "; e_in = " << e_in << "; e_out = " << e_out << 
  		"\nionic_strenght = " << ionic_strength << "; decay = " << decay << std::endl;
  std::cout << "Linear solver = " << linear_solver_name << std::endl;
  std::cout << "Linear solver options = " << linear_solver_options << "\n" << std::endl;
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
poisson_boltzmann::is_in (const NS::Atom& i,
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

  retval =           (i.pos[0] > l - tol) && (i.pos[0] <= r  - tol);
  retval = retval && (i.pos[1] > f - tol) && (i.pos[1] <= bk - tol);
  retval = retval && (i.pos[2] > b - tol) && (i.pos[2] <= t  - tol);

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
            double min = 100.0 * this->atoms.size (); 
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
                  for (const NS::Atom& i : atoms) 
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
            double min = 100.0 * this->atoms.size (); 
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

                for (const NS::Atom& i : atoms) 
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
  this->rho_fixed.assign (this->tmsh.num_local_quadrants (), 0.0);  
  
  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
       quadrant != this->tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (const NS::Atom& i : atoms) //(const ion& i : ions)
        if (is_in (i, quadrant)) //Se un atomo è nel quadrante che sto considerando:
          {
            double volume = (quadrant->p(0, 7) - quadrant->p(0, 0)) *
              (quadrant->p(1, 7) - quadrant->p(1, 0)) *
              (quadrant->p(2, 7) - quadrant->p(2, 0)); //calcolo il volume del quadrante 
            //this->rho_fixed[quadrant->get_forest_quad_idx ()] = i.charge / volume; 
            //assegno carica_corrispondente/volume al quadrante come valore di rho corrispondente
            //rendo adimensionale:
            this->rho_fixed[quadrant->get_forest_quad_idx ()] = -(i.charge / volume)*4.0*pi;  //*Angs*Angs*Angs/e;
            //con il meno se è a dx 
            break;
          }

      int num_int_nodes = 0;
      int num_hanging = 0;
      for (int ii = 0; ii < 8; ++ii)
        {
          if (! quadrant->is_hanging (ii)) //se il quadrante corrente NON è stato "appeso"
            {
              if (this->levelsetfun (quadrant->p (0, ii),
                                     quadrant->p (1, ii),
                                     quadrant->p (2, ii)) > 1.0)
                //restituisce una distanza per capire se sono dentro o fuori dalla molecola. Se > 1.0 -> Dentro
                ++num_int_nodes; //aumento numero dei nodi dentro la molecola 
            }
          else
            ++num_hanging; //aumento numero hanging 
        }
      if (num_int_nodes == 0)  //se non ho nodi dentro
        this->marker[quadrant->get_forest_quad_idx ()] = 1.0; //assegno al marker del quadrante 1.0 -> fuori
      else if (num_int_nodes < (8 - num_hanging)) //se numero di nodi interni è meno di quelli fuori 
        this->marker[quadrant->get_forest_quad_idx ()] = 1.0/2.0; //assegno a quel quad 1/2 
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
  // diffusion
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/(e*e);   //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/(e*e); //adim e_out
  epsilon.assign (tmsh.num_local_quadrants (), eps_in); //e_in
  
  //reaction 
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/(e_0*e_out*kb*T); 
  
  for (auto epsp = epsilon.begin (), mp = marker.begin ();
       epsp != epsilon.end () || mp != marker.end ();
       ++epsp, ++mp) //scorro vettore delle epsilon e dei marker
    if ((*mp) == 0.0) //se marker vale zero -> sono dentro -> e_in
      (*epsp) = eps_in; //e_in
    else //se marker è diverso da zero -> e_out 
      (*epsp) = eps_out; //e_out

  reaction.assign (tmsh.num_local_quadrants (), 0.0);
  for (auto rp = reaction.begin (), mp = marker.begin ();
       rp != reaction.end () || mp != marker.end ();
       ++rp, ++mp)
    if ((*mp) != 0.0) //se marker è diverso da zero -> sono fuori (o il quad è più fuori che dentro) -> ho k^2*eps_out
      (*rp) = -eps_out*k2; //dentro è sempre 0.0, con il meno se è a sx

  tmsh.octbin_export_quadrant ("epsilon_0", epsilon);
  tmsh.octbin_export_quadrant ("rho_0", rho_fixed);
  tmsh.octbin_export_quadrant ("reaction_0", reaction);

  distributed_sparse_matrix A; //matrice A per diffusione 
  A.set_ranges (tmsh.num_owned_nodes ());
  
  A.resize (tmsh.num_global_nodes ());
  
  distributed_vector  rhs (tmsh.num_global_nodes ()); //vettore del rhs

  distributed_vector  psi (tmsh.num_global_nodes ());
  psi.get_owned_data ().assign (psi.get_owned_data ().size (), 0.0); //creo psi nulla (no trasporto)
  
  bim3a_solution_with_ghosts (tmsh, psi);
  
  distributed_vector ones (tmsh.num_global_nodes ());
  ones.get_owned_data ().assign (ones.get_owned_data ().size (), 1.0); //vettore di 1.0
  
  bim3a_solution_with_ghosts (tmsh, ones, replace_op);
  
  bim3a_advection_diffusion (tmsh, epsilon, psi, A); //diffusion con i suoi coeff epsilon (e_in dentro, e_out fuori)
  bim3a_reaction (tmsh, reaction, ones, A); //reaction (reazione è 0 dentro, fuori eps*k2)
  
  A.assemble (); //Assemblo la matrice con diffusione e reazione 

  bim3a_rhs (tmsh, rho_fixed, ones, rhs); //rhs (vale rho_fixed o zero)
  
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
  
  /*
  std::ofstream fout ((std::string ("Sol_mumps.txt")).c_str ());
  fout << std::endl;

  for (unsigned int k = 0; k < phi.size (); ++k)
  	fout << "mumps: " << phi[k] << std::endl;
  
  fout.close ();
  
  
  //lis:
  std::cout << "\nStarting lis solution" << std::endl; 
  
  sparse_matrix A_lis;
  A_lis.resize (tmsh.num_global_nodes ());
  std::vector<double> rhs_lis(tmsh.num_global_nodes (), 0.0);
  std::vector<double> xa;
  
  std::cout << "\nSolver" << std::endl; 
  lis lis_solver;
  lis_solver.set_iterative_method("Conjugate Gradient"); //metodo del gradiente coniugato
  //lis_solver.set_tolerance (1e-14); //tolleranza 
  //lis_solver.set_max_iterations (500); //number of iterations
  lis_solver.set_preconditioner ("jacobi"); //non credo ci sia Cholesky, c'è lu o ilu o jacobi
  
  int base = lis_solver.get_index_base ();
  
  std::string prec;
  lis_solver.get_preconditioner (prec);
  std::cout << "Using solver of type "
            << lis_solver.solver_type ()
            << " named "
            << lis_solver.solver_name ()
            << " with preconditioner "
            << prec
            << std::endl;  
            
  std::vector<double> psi_lis(tmsh.num_local_quadrants (), 0.0); 
  std::vector<double> ones_lis(tmsh.num_local_quadrants (), 1.0);
 
  bim3a_structure (tmsh, A_lis);
  
  bim3a_advection_diffusion (tmsh, epsilon, psi_lis, A_lis);
  
  bim3a_reaction (tmsh, reaction, ones_lis, A_lis);
  
  bim3a_rhs (tmsh, rho_fixed, ones_lis, rhs_lis);
  
  A_lis.aij (xa, irow, jcol, base);
  
  linear_solver::matrix_format_t mf = linear_solver::matrix_format_t::aij; //csr o aij??
  lis_solver.set_lhs_structure (A_lis.rows (), irow, jcol, mf);
  lis_solver.set_lhs_data (xa);
  lis_solver.set_rhs (rhs_lis);
  
  std::cout << "\nlis_solver.analyze () = "
            << lis_solver.analyze ()
            << std::endl;
  std::cout << "lis_solver.factorize () = "
            << lis_solver.factorize ()
            << std::endl;
  std::cout << "lis_solver.solve () = "
            << lis_solver.solve ()
            << std::endl;
  
  std::cout << "\nResult of PBE with lis \nwill be written in "
            << "Compare_solution.txt"
            << std::endl;

  std::ofstream fout ((std::string ("Compare_solution.txt")).c_str ());
  fout << std::endl;

  for (unsigned int k = 0; k < rhs_lis.size (); ++k)
  	fout << "mumps: " << rhs[k] 
  	     << "	lis: " << rhs_lis[k] 
  	     << "	diff: " << rhs[k] - rhs_lis[k] << std::endl;
  
  fout.close ();

  lis_solver.cleanup ();
  */
  
  // CON DISTRIBUTED STRUCTURES:
  //distributed_sparse_matrix A_lis;
  //A_lis.set_ranges (tmsh.num_owned_nodes ());
  //A_lis.resize (tmsh.num_global_nodes ());
  //distributed_vector  rhs_lis (tmsh.num_global_nodes ()); //vettore del rhs

}
