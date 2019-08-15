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


struct
poisson_boltzmann
{

  struct
  ion
  {
    std::array<double, 3> center;
    double radius;
    double charge;
    
    ion (const std::array<double, 3>& c,
         const double& r)
      : center (c), radius (r) {};
  };

  static constexpr p4est_topidx_t simple_conn_num_vertices = 8;
  static constexpr p4est_topidx_t simple_conn_num_trees = 1;
  std::array <double, simple_conn_num_vertices*3>  simple_conn_p;
  std::array <p4est_topidx_t, simple_conn_num_trees*9>  simple_conn_t;
  
  std::vector<ion> ions;

  double ll;
  double rr;
  
  int maxlevel;
  int minlevel;
  double decay;
  double e_in, e_out;
  
  MPI_Comm mpicomm;  
  tmesh_3d tmsh;
  
  std::string csvfilename;
  std::string p4estfilename;
  std::string lsfilename;
  std::string markerfilename;

  std::vector<double> marker;
  std::vector<double> epsilon;
  std::vector<double> rho_fixed;
  std::vector<double> reaction;
  
  poisson_boltzmann (int maxlevel_ = 8, int minlevel_ = 3,
                     double decay_ = -1.5, double e_in_ = 4.0,
                     double e_out_ = 80.0,
                     MPI_Comm mpicomm_ = MPI_COMM_WORLD)
    : maxlevel(maxlevel_),
      minlevel(minlevel_),
      decay(decay_),
      e_in(e_in_),
      e_out(e_out_),
      mpicomm(mpicomm_),
      tmsh(mpicomm)
  {  };
    
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
  is_in (const ion& i, tmesh_3d::quadrant_iterator q);
  
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
