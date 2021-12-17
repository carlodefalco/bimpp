#ifndef HAVE_PB_CLASS_H
#define HAVE_PB_CLASS_H

#include <bim_timing.h>
#include <tmesh_3d.h>

#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }
#include <cmath>

const double p4esttol = 1 / std::pow (2, P8EST_QMAXLEVEL);

#include <array>
#include <string>
#include <vector>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//#include "pqr_parser.cpp"
#include "nanoshaper.h"

// Problem parameters
constexpr double e_0 = 8.85418781762e-12;	//Dielectric void const [F/m]
constexpr double kb = 1.380649e-23;		//Boltzmann constant [J/K]
constexpr double T = 273.15 + 25;		//Temperature [K]
constexpr double e = 1.602176634e-19;  	//Charge of an electron [C]
constexpr double N_av = 6.022e23;      	//Avogadro Number [mol^-1]
constexpr double Angs = 1e-10;         	//Angstrom [m]
constexpr double pi = 3.14159265358979323846; 

struct
poisson_boltzmann
{

  static constexpr p4est_topidx_t simple_conn_num_vertices = 8;
  static constexpr p4est_topidx_t simple_conn_num_trees = 1;
  std::array <double, simple_conn_num_vertices*3>  simple_conn_p;
  std::array <p4est_topidx_t, simple_conn_num_trees*9>  simple_conn_t;
  
  std::vector<NS::Atom> atoms;

  //Cubic mesh:
  double ll; //min value between all the coordinates 
  double rr; //max value between all the coordinates
  
  //Stretched mesh:
  double l_c[3]; //min x, y, z value
  double r_c[3]; //max x, y, z value
  
  //mesh:
  int maxlevel;
  int minlevel;
  int mesh_shape;
  
  //model:
  int linearized;
  double decay;
  double e_in, e_out, ionic_strength; //[M]
  
  //algorithm:
  std::string linear_solver_name;
  std::string linear_solver_options;
  std::string linear_solver_preconditioner;
  std::string linear_solver_precond_opts;
  std::string linear_solver_tol;

  MPI_Comm mpicomm;
  tmesh_3d tmsh;

  std::string optionsfilename;
  std::string pqrfilename;
  std::string p4estfilename;
  std::string lsfilename; 
  std::string markerfilename;

  std::vector<double> marker; 
  std::vector<double> epsilon; 
  std::vector<double> rho_fixed; 
  std::vector<double> reaction; 

  poisson_boltzmann (int maxlevel_ = 4, int minlevel_ = 3, int mesh_shape_ = 1,
                     int linearized_ = 1, double decay_ = -1.5,
                     double e_in_ = 2.0, double e_out_ = 80.0, double ionic_strength_ = 0.145,
                     std::string linear_solver_name_ = "mumps", std::string linear_solver_options_ = "",
                     MPI_Comm mpicomm_ = MPI_COMM_WORLD)
    : maxlevel(maxlevel_),
      minlevel(minlevel_),
      mesh_shape(mesh_shape_),
      linearized(linearized_),
      decay(decay_),
      e_in(e_in_),
      e_out(e_out_),
      ionic_strength(ionic_strength_),
      linear_solver_name(linear_solver_name_),
      linear_solver_options(linear_solver_options_),
      mpicomm(mpicomm_),
      tmsh(mpicomm)
  {  };

  double
  levelsetfun (double x, double y, double z);

  static int
  uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
  { return 1; }
  
  void
  create_cubic_mesh ();
  
  void
  create_mesh ();
  
  int
  parse_options (int argc, char **argv);
  
  void 
  print_options ();

  void
  init_tmesh ();

  bool
  is_in (const NS::Atom& i, tmesh_3d::quadrant_iterator q); 

  void
  refine_surface ();

  void
  create_markers ();

  void
  export_ls_tmesh ();

  void
  export_marked_tmesh ();

  void
  export_p4est ();

  void
  mumps_compute_electric_potential ();
  
  void
  lis_compute_electric_potential ();

};

#endif
