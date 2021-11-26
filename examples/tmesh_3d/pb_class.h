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
/*
  struct
  ion
  {
    std::array<double, 3> center;
    double radius;
    double charge;

    ion (const std::array<double, 3>& c,
         const double& r)
      : center (c), radius (r) {};
  };*/

  static constexpr p4est_topidx_t simple_conn_num_vertices = 8;
  static constexpr p4est_topidx_t simple_conn_num_trees = 1;
  std::array <double, simple_conn_num_vertices*3>  simple_conn_p;
  std::array <p4est_topidx_t, simple_conn_num_trees*9>  simple_conn_t;

  //std::vector<ion> ions;
  
  std::vector<NS::Atom> atoms;

  double ll;
  double rr;

  int maxlevel;
  int minlevel;
  double decay;
  double e_in, e_out, k2;

  MPI_Comm mpicomm;
  tmesh_3d tmsh;

  std::string csvfilename;
  std::string p4estfilename;
  std::string lsfilename;
  std::string markerfilename;

  std::vector<double> marker; //vettore che mi dice se sono dentro o fuori dalla molecola
  std::vector<double> epsilon; //Vettore che vale e_in se dentro molecola, e_out altrimenti
  std::vector<double> rho_fixed; //vettore delle cariche fisse
  std::vector<double> reaction; //vettore del termine di reazione: eps(r)*k^2 (k=A^2/lambda^2)
  //k2 è nullo dentro la molecola e nello stern layer 

  poisson_boltzmann (int maxlevel_ = 4, int minlevel_ = 3, //maxlevel_ = 8, minlevel_ = 3
                     double decay_ = -1.5, double e_in_ = 2.0, //e_in_ = 4.0
                     double e_out_ = 80.0, double k2_ = 1.0,
                     MPI_Comm mpicomm_ = MPI_COMM_WORLD)
    : maxlevel(maxlevel_),
      minlevel(minlevel_),
      decay(decay_),
      e_in(e_in_),
      e_out(e_out_),
      k2(k2_),
      mpicomm(mpicomm_),
      tmsh(mpicomm)
  {  };
  //minlevel = numero di raffinamenti uniformi (è nella funzione init mesh)
  //maxlevel - minlevel = numero di raffinamenti adattivi (in refine surface) ??

  double
  levelsetfun (double x, double y, double z);

  static int
  uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
  { return 1; }

  void
  read_csv ( );
  
  void
  parse_options (int argc, char **argv);

  void
  init_tmesh ();

  bool
  is_in (const NS::Atom& i, tmesh_3d::quadrant_iterator q); //(const ion& i, tmesh_3d::quadrant_iterator q);

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
  compute_electric_potential ();

};

#endif
