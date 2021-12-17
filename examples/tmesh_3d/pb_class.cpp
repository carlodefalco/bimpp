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
            this->rr = a1.pos[kk]; 
          else if (a1.pos[kk] < this->ll)
            this->ll = a1.pos[kk]; 
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
            this->l_c[kk] = a1.pos[kk]; 
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
  optionsfilename = g ("potfile", "../options.pot");
  pqrfilename = g ("pqrfile", "1CCM.pqr");
  
  //Check that the input files exist 
  std::ifstream optionsfile(optionsfilename);
  if(!optionsfile)
  {
  	std::cerr << "Cannot find the options file" << std::endl;
  	return 1;
  }
  
  std::ifstream pqrfile(pqrfilename);
  if(!pqrfile)
  {
  	std::cerr << "Cannot find the pqr file" << std::endl;
  	return 1;
  }
  
  //Read the options from the file
  GetPot g2 (optionsfilename.c_str ());
  
  const std::string mesh_options = "mesh/";
  maxlevel = g2 ((mesh_options + "maxlevel").c_str (),  6);
  minlevel = g2 ((mesh_options + "minlevel").c_str (),  4);
  mesh_shape = g2 ((mesh_options + "mesh_shape").c_str (),  1);

  const std::string model_options = "model/";
  linearized = g2 ((model_options + "linear_solver").c_str (),  1);
  ionic_strength = g2 ((model_options + "ionic_strength").c_str (),  0.145);
  e_in = g2 ((model_options + "molecular_dielectric_constant").c_str (),  2.);
  e_out = g2 ((model_options + "solvent_dielectric_constant").c_str (),  78.54);
  decay    = g2 ((model_options + "decay").c_str (), -1.5);
  
  const std::string alg_options = "algorithm/";
  linear_solver_name = g2 ((alg_options + "linear_solver").c_str (),  "mumps");
  linear_solver_options = g2 ((alg_options + "solver_options").c_str (),  "cg");
  linear_solver_preconditioner = g2 ((alg_options + "preconditioner").c_str (),  "ilu");
  linear_solver_precond_opts = g2 ((alg_options + "preconditioner_options").c_str (),  "0");
  linear_solver_tol = g2 ((alg_options + "tol").c_str (),  "1.e-12");

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
  std::cout << "minlevel = " << minlevel <<  "\nmaxlevel = " << maxlevel << std::endl;
  if(mesh_shape == 1)
  	std::cout << "Mesh shape = stretched" << std::endl;
  else
  	std::cout << "Mesh shape = cubic" << std::endl;
  std::cout << "Linearized model = " << linearized << "\ne_in = " << e_in << "\ne_out = " << e_out << 
  		"\nionic_strenght = " << ionic_strength << "\ndecay = " << decay << std::endl;
  std::cout << "Linear solver = " << linear_solver_name << std::endl;
  std::cout << "Linear solver options = " << linear_solver_options << std::endl;
  std::cout << "Preconditioner = " << linear_solver_preconditioner << std::endl;
  if (linear_solver_preconditioner == "ilu")
     std::cout << "ilu fill level selected = " <<  linear_solver_precond_opts << std::endl;
  std::cout << "Choosen tolerance = " << linear_solver_tol << "\n" << std::endl;
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
      for (const NS::Atom& i : atoms) 
        if (is_in (i, quadrant)) 
          {
            double volume = (quadrant->p(0, 7) - quadrant->p(0, 0)) *
              (quadrant->p(1, 7) - quadrant->p(1, 0)) *
              (quadrant->p(2, 7) - quadrant->p(2, 0)); //volume
            this->rho_fixed[quadrant->get_forest_quad_idx ()] = -(i.charge / volume)*4.0*pi;
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
poisson_boltzmann::mumps_compute_electric_potential ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  if (rank == 0)
     std::cout << "\nStarting MUMPS solution" << std::endl;
     
  // diffusion
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/(e*e);   //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/(e*e); //adim e_out
  epsilon.assign (tmsh.num_local_quadrants (), eps_in); //e_in
  
  //reaction 
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/(e_0*e_out*kb*T); 
  
  for (auto epsp = epsilon.begin (), mp = marker.begin ();
       epsp != epsilon.end () || mp != marker.end ();
       ++epsp, ++mp) 
    if ((*mp) == 0.0)
      (*epsp) = eps_in; 
    else  
      (*epsp) = eps_out; 

  reaction.assign (tmsh.num_local_quadrants (), 0.0);
  for (auto rp = reaction.begin (), mp = marker.begin ();
       rp != reaction.end () || mp != marker.end ();
       ++rp, ++mp)
    if ((*mp) != 0.0) 
      (*rp) = -eps_out*k2; 

  tmsh.octbin_export_quadrant ("epsilon_0", epsilon);
  tmsh.octbin_export_quadrant ("rho_0", rho_fixed);
  tmsh.octbin_export_quadrant ("reaction_0", reaction);

  distributed_sparse_matrix A;  
  A.set_ranges (tmsh.num_owned_nodes ());
  
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
  
  A.assemble (); 
  rhs.assemble();

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


void
poisson_boltzmann::lis_compute_electric_potential ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  if (rank == 0)
     std::cout << "\nStarting LIS solution" << std::endl;
     
  // diffusion
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/(e*e);   //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/(e*e); //adim e_out
  epsilon.assign (tmsh.num_local_quadrants (), eps_in); //e_in
  
  //reaction 
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/(e_0*e_out*kb*T); 
  
  for (auto epsp = epsilon.begin (), mp = marker.begin ();
       epsp != epsilon.end () || mp != marker.end ();
       ++epsp, ++mp) 
    if ((*mp) == 0.0)
      (*epsp) = eps_in; 
    else  
      (*epsp) = eps_out; 

  reaction.assign (tmsh.num_local_quadrants (), 0.0);
  for (auto rp = reaction.begin (), mp = marker.begin ();
       rp != reaction.end () || mp != marker.end ();
       ++rp, ++mp)
    if ((*mp) != 0.0) 
      (*rp) = -eps_out*k2; 

  tmsh.octbin_export_quadrant ("epsilon_0", epsilon);
  tmsh.octbin_export_quadrant ("rho_0", rho_fixed);
  tmsh.octbin_export_quadrant ("reaction_0", reaction);
  
  distributed_sparse_matrix A; 
  A.set_ranges (tmsh.num_owned_nodes ());
  
  //A.resize (tmsh.num_global_nodes ()); //old
  
  //distributed_vector  rhs (tmsh.num_global_nodes ()); //old
  distributed_vector  rhs (tmsh.num_owned_nodes (), mpicomm); //new

  distributed_vector  psi (tmsh.num_global_nodes ());
  psi.get_owned_data ().assign (psi.get_owned_data ().size (), 0.0);

  bim3a_solution_with_ghosts (tmsh, psi);
  
  distributed_vector ones (tmsh.num_global_nodes ()); 
  ones.get_owned_data ().assign (ones.get_owned_data ().size (), 1.0); 
  
  bim3a_solution_with_ghosts (tmsh, ones, replace_op);
  
  bim3a_advection_diffusion (tmsh, epsilon, psi, A);
  bim3a_reaction (tmsh, reaction, ones, A); 

  A.assemble (); 
  rhs.assemble();

  bim3a_rhs (tmsh, rho_fixed, ones, rhs);
  
  tmsh.octbin_export ("rhs_0", rhs);  
  
  //CSR  
  std::vector<double> vals;
  std::vector<int> irow, jcol;

  A.csr(vals, jcol, irow);

  // lis RHS
  LIS_INT i, is, ie, n_rhs, ln; 
  LIS_VECTOR rhs_lis;
  n_rhs = tmsh.num_global_nodes();
  //ln = tmsh.num_local_nodes();
  ln = rhs.get_owned_data().size(); //equivalent to: tmsh.num_owned_nodes();
  
  lis_vector_create(mpicomm, &rhs_lis);
  lis_vector_set_size(rhs_lis, ln, 0);
  //lis_vector_set_size(rhs_lis, 0, n_rhs); 
  lis_vector_get_range(rhs_lis, &is, &ie);
  
  lis_vector_set_values2(LIS_INS_VALUE, is, ln, &(rhs.get_owned_data()[0]), rhs_lis); //pass values to rhs_lis 
  //lis_vector_print(rhs_lis);
  
  // lis PHI
  LIS_VECTOR phi_lis;
  
  lis_vector_create(mpicomm, &phi_lis);
  lis_vector_set_size(phi_lis, ln, 0);
  //lis_vector_set_size(phi_lis, 0, n_rhs);
  lis_vector_get_range(phi_lis, &is, &ie); //is, ie assigned again
  
  // lis MATRIX
  LIS_INT n, nnz; //n: matrix dim ; nnz: numb of non zero elems 
  LIS_INT *index; //array of integer containing the col index of non zero elems
  LIS_INT *ptr; //array of integer with starting points of rows 
  LIS_SCALAR *value; //array of double stores non-zero elements of matrix A along the row
  LIS_MATRIX A_lis; //array of integer containing the col index of non zero elems 
  
  nnz = A.owned_nnz();
  n = tmsh.num_owned_nodes ();
  
  //Some prints:
  /*
  std::cout << "n_global = " << tmsh.num_global_nodes() << std::endl;
  std::cout << "n_local = " << tmsh.num_local_nodes() << std::endl;
  std::cout << "nnz = " << nnz << std::endl;
  std::cout << "n = " << n << std::endl;
  std::cout << "vals size: " << vals.size() << std::endl;
  std::cout << "jcol size: " << jcol.size() << std::endl;
  std::cout << "irow size: " << irow.size() << std::endl;
  */
    
  //ptr = (LIS_INT *)malloc( (n+1)*sizeof(LIS_INT) );
  //index = (LIS_INT *)malloc( nnz*sizeof(LIS_INT) );
  //value = (LIS_SCALAR *)malloc( nnz*sizeof(LIS_SCALAR) );
  
  lis_matrix_create(mpicomm, &A_lis);
  
  lis_matrix_set_size(A_lis, n, 0);

  ptr = &irow[0];
  index = &jcol[0];
  value = &vals[0];
  
  lis_matrix_set_csr(nnz, ptr, index, value, A_lis);
  
  lis_matrix_assemble(A_lis);
  
  //Solve linear system
  LIS_SOLVER solver;
  
  lis_solver_create(&solver);
  
  std::string opts = "-i " + linear_solver_options + " -p " + linear_solver_preconditioner;
  if (linear_solver_preconditioner == "ilu" && linear_solver_precond_opts != "")
     opts += " -ilu_fill [" + linear_solver_precond_opts + "]";
     
  opts += " -tol " + linear_solver_tol;
  
  //std::cout << "opts : " << opts << std::endl;
  lis_solver_set_option(&opts[0], solver);
  
  //LIS_INT s, p;
  //lis_solver_get_solver(solver, &s);
  //lis_solver_get_precon(solver, &p);
  //std::cout << "Solver = " << s << " ; precon = " << p << std::endl;
  
  lis_solve(A_lis, rhs_lis, phi_lis, solver);

  distributed_vector phi (tmsh.num_owned_nodes (), mpicomm); 

  lis_vector_get_values(phi_lis, is, ln, &(phi.get_owned_data ()[0]));

  bim3a_solution_with_ghosts (tmsh, phi);

  tmsh.octbin_export ("phi_0", phi); 
  
}
