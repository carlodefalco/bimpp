/*
  Copyright (C) 2020 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include <octave_file_io.h>
 
#include <bim_distributed_vector.h>
#include <bim_timing.h>
#include <tmesh.h>
#include <quad_operators.h>

#include "Taylor_Galerkin.h"

static constexpr char LOADFILENAME_1[255] = "orography.octbin.gz";
static constexpr char LOADFILENAME_2[255] = "basin.octbin.gz";
static constexpr char SAVEFILENAME_1[255] = "orography_tmsh";
static constexpr char VARNAME_1[255] = "dem";
static constexpr char VARNAME_2[255] = "basin";

// properties of the input dem
static constexpr double res = 5; // it is also the minimum resolution of the bim element
static constexpr double Nx = 1998;
static constexpr double Ny = 1829;


static constexpr double L = 2;//res*(Nx-1);
static constexpr double H = 2;//res*(Ny-1);
static std::vector<double>   dem;
static std::vector<double>   basin_mask;
static constexpr int NUM_REFINEMENTS = 6;
static constexpr int NUM_TREFINEMENTS = 1; // 10


static constexpr double SAVEDT  = 1e-3;
static constexpr double DELTAT =  1e-3;
static constexpr double REDCDT =  1;
static constexpr double T      =  5.;





// Connectivity of local element
constexpr p4est_topidx_t simple_conn_num_vertices = 4;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*2] =
  {0,  0,
   0,  H,
   L,  0,
   L,  H};

const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
  {  1,    3,    4,    2,    1 };

// Refinement rule
static int
uniform_refinement (tmesh::quadrant_iterator q)
{ return NUM_REFINEMENTS; }

static int
hanging_refinement (tmesh::quadrant_iterator quadrant)
{
  double x_minus, x_plus, y_minus, y_plus;
  x_minus = quadrant->p (0, 0);
  x_plus  = quadrant->p (0, 1);
  y_minus = quadrant->p (1, 0);
  y_plus  = quadrant->p (1, 2);
  
  double x_center, y_center;
  x_center = quadrant->centroid (0);
  y_center = quadrant->centroid (1);
  
  const auto marker = x_center<3./4.*L && x_center>L/4. && y_center<3./4.*H && y_center>H/4.;
//  const auto marker = (x_minus+x_plus)/2<L/2;//(x_minus+x_plus)/2<L/2 && (y_minus+y_plus)/2>H/2 ? 1 : 0;
  return marker; }


static int
raster_2_vector(const double& i_x,
                const double& i_y)
{
  static int ii;
  ii = i_y + Ny*i_x;
  
  return(ii);
}

static std::array<int,3>
global_coord_2_raster(const double& x,
                      const double& y)
{
  static double i_x;
  static double i_y;
  static int ii;
  
  
  // nearest neighbor
  i_x = std::round(x / res);
  i_y = - std::round(y / res) + (Ny-1);
  
  ii = raster_2_vector(i_x,i_y);
  
  return(std::array<int,3>{{ ii,int(i_x),int(i_y) }});
}



using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector
using Q0  = std::vector<double>;         // Typedef for local q_0 vector

//double h0_fun (const double& xx, const double& yy)  { return std::max (0., (8. - std::sin (M_PI * xx / 2. / 400.) - dem[global_coord_2_raster(xx,yy)[0]])); }
double h0_fun (const double& xx, const double& yy, const double& L, const double& H)  {
  
  
//  return ( 1.+1.*std::exp(-0.5*( std::pow(xx-L/2.,2.)+std::pow(yy-H/2.,2.) )/std::pow(0.2*L/2.,2.) ) );
//  return ( 0.+1.*std::exp(-0.5*( std::pow(xx-L/2.,2.)+std::pow(yy-H/2.,2.) )/std::pow(0.2*L/2.,2.) ) );

  
  
  if (xx>L*1./4. && xx<L*3./4. && yy >H*1./4. && yy <H*3./4.) //(xx>L*3./10. && xx<7./10.*L)
  {
    return 1.;
  }
  
//  if (xx>L/4 && xx<3/4*L && yy>H/4 && yy <3/4*H)
//  {
//    return 2;
//  }
  return 0.;
}
double Ux0_fun (double xx, double yy) { return 0.; }
double Uy0_fun (double xx, double yy) { return 0.; }

// Assemble vector from mesh.
// FIXME  the following two functions are copied over from
// "quad_operators.cpp" as they were not exported in an header,
// should find better way to avoid code duplication
void
assemble_vector (tmesh::quadrant_iterator& quadrant,
                 const std::array<double, 4>& locrhs,
                 Q1& rhs,
                 const ordering& ord = default_ord)
{
  
  std::vector<unsigned int> rows;
  rows.reserve (2);
  int i, r;
  
  for (i = 0; i < 4; ++i)
  {
    rows.clear ();
    
    if (! quadrant->is_hanging (i))
      rows.push_back (quadrant->gt (i));
    else
    {
      rows.push_back (quadrant->gparent (0, i));
      rows.push_back (quadrant->gparent (1, i));
    }
    
    for (r = 0; r < rows.size (); ++r)
    rhs[ord (rows[r])] +=
    locrhs[i] / rows.size ();
  }
}





// Re-Define tic and toc to add an MPI_Barrier
#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }


int
main (int argc, char **argv)
{
  // Management of solutions ordering
  ordering ordh  = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 0> (gt); };
  ordering ordUx = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 1> (gt); };
  ordering ordUy = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<3, 2> (gt); };

  
  
  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  
  
  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  
  MPI_Barrier (MPI_COMM_WORLD);
  
  
  TIC ();
  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);
  
  
//  tmsh.set_refine_marker (hanging_refinement);
//  tmsh.refine (recursive, 1);
  TOC ("Uniform refinement");
  
  
  
  // ln_nodes sono i dof non gli hanging node!! (sono esclusi dal calcolo)
  tmesh::idx_t gn_nodes    = tmsh.num_global_nodes (); // Return total number of nodes owned by all process
  tmesh::idx_t ln_nodes    = tmsh.num_owned_nodes (); // Return number of nodes owned by local process
  tmesh::idx_t ln_elements = tmsh.num_local_quadrants ();  // Return number of quadrants owned by local process across all trees
  tmesh::idx_t gn_elements = tmsh.num_global_quadrants (); // Return number of quadrants owned by all processes across all trees
  
  
  /// Allocate initial data container
  Q1 sol  (ln_nodes * 3);
  Q1 incr (ln_nodes * 3);
  sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);
  incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);
  
  Q1 mass (ln_nodes * 3);
  bim2a_mass_vector (tmsh, mass, ordh);
  bim2a_mass_vector (tmsh, mass, ordUx);
  bim2a_mass_vector (tmsh, mass, ordUy);
  mass.assemble ();
  
  Q0 sol_onehalf(gn_elements * 3);
  sol_onehalf.assign(sol_onehalf.size(), 0.0);
  
  
  // .assemble serve probabilmente in un ambiente parallelo per creare la matrice di massa globale
  
//   gioca un po' con questo per capire la matrice di massa lumpata e come integrare il termine di flusso
//  for (const auto& it : mass.get_owned_data ())
//  {
//    std::cout << it << std::endl;
//  }
//
//  std::cout << "STOP!" << " " << ln_nodes << " " << ln_elements << " " << mass.get_owned_data ().size() <<std::endl;
//  return 0;
  
  
  Q1 Z    (ln_nodes);
  Z.get_owned_data ().assign (Z.get_owned_data ().size (), 0.0);
  

  TIC();
//  octave_io_mode m_in = gz_read_mode, m_out = gz_read_mode;
//  octave_value v;
//
//  octave_io_open (LOADFILENAME_1, m_in, &m_out);
//  octave_load (VARNAME_1, v);
//  Matrix M = v.matrix_value ();
//  dem.resize (M.numel ());
//  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), dem.begin ());
//
//  octave_io_open (LOADFILENAME_2, m_in, &m_out);
//  octave_load (VARNAME_2, v);
//  M = v.matrix_value ();
//  basin_mask.resize (M.numel ());
//  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), basin_mask.begin ());
  
  TOC("Load data matrix");
  
  
  // Buffer for export filename
  char filename[255]="";
  

  // Initialize initial datas
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
  {
    
    for (int ii = 0; ii < 4; ++ii)
    {
      if (! quadrant->is_hanging (ii)){
        double xx=quadrant->p(0,ii);
        double yy=quadrant->p(1,ii);
        
        sol [ordh (quadrant->gt (ii))] = h0_fun (xx, yy, L, H);
        sol [ordUx(quadrant->gt (ii))] = Ux0_fun (xx, yy);
        sol [ordUy(quadrant->gt (ii))] = Uy0_fun (xx, yy);
        
        Z[quadrant->gt (ii)] = 0.;//dem[global_coord_2_raster(xx,yy)[0]]*0;
        
//        if (h0_fun (xx, yy) !=0) std::cout << h0_fun (xx, yy) << std::endl;
      }
      
      else
      {
        // touch parent nodes to set up distributed vector structure
        
        sol [ordh (quadrant->gparent(0,ii))] += 0.;
        sol [ordh (quadrant->gparent(1,ii))] += 0.;
        sol [ordUx(quadrant->gparent(0,ii))] += 0.;
        sol [ordUx(quadrant->gparent(1,ii))] += 0.;
        sol [ordUy(quadrant->gparent(0,ii))] += 0.;
        sol [ordUy(quadrant->gparent(1,ii))] += 0.;
        
        Z   [quadrant->gparent(0,ii)] += 0.;
        Z   [quadrant->gparent(1,ii)] += 0.;
      }
    }
  }
  // bim2a_solution_with_ghosts in quad_operators.cpp
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, sol, replace_op, ordUy);
  
  bim2a_solution_with_ghosts (tmsh, Z, replace_op);
  
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordh,  false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUx, false);
  bim2a_solution_with_ghosts (tmsh, incr, replace_op, ordUy);
  TOC ("compute initial condition");
  
  // create results dir
//  if (rank == 0)
//  {
//    std::string bashCommand = std::string ( "rm -r results" );
//    std::system ( bashCommand.c_str() );
//
//    bashCommand = std::string ( "mkdir -p results" );
//    std::system ( bashCommand.c_str() );
//  }

  
  
  // Save initial conditions
  sprintf(filename, "results/swe_h_%4.4d", 0);
  tmsh.octbin_export (filename, sol, ordh);
  sprintf(filename, "results/swe_Ux_%4.4d", 0);
  tmsh.octbin_export (filename, sol, ordUx);
  sprintf(filename, "results/swe_Uy_%4.4d", 0);
  tmsh.octbin_export (filename, sol, ordUy);
  sprintf(filename, "results/swe_Z_%4.4d", 0);
  tmsh.octbin_export (filename, Z);
  

  std::vector<double> full_time_vector;
  full_time_vector.reserve (static_cast<int> (T/DELTAT));
  std::vector<double> save_time_vector;
  save_time_vector.reserve (static_cast<int> (T/SAVEDT));
  
  
  
  TG2_scheme stp(sol, sol_onehalf, ordh, ordUx, ordUy, Z, DELTAT);
  
  
  
  // Time loop
  double time = 0.0;
  double deltat = DELTAT;
  if(rank==0) {
    full_time_vector.push_back (0.0);
    save_time_vector.push_back (0.0);
  }
  int count = 0;
  
  double savecount = 0.0;
  
  if (rank == 0)
  {
    std::cout << "start loop" << std::endl;
  }
  
  while (time <= T)
  {
    // Reset increment
    TIC();
    incr.get_owned_data ().assign (incr.get_owned_data ().size (), 0.0);
    incr.assemble (replace_op);
    TOC("Reset");
    TIC();
    
  
    // compute time step,
    stp.set_dt (DELTAT);
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.compute_dt(quadrant);
    }
    deltat = REDCDT * stp.dt;
    

    
    MPI_Allreduce (MPI_IN_PLACE, static_cast<void*> (&deltat), 1, MPI_INT, MPI_MIN, tmsh.comm);
    time += deltat;
    savecount += deltat;
    MPI_Bcast (static_cast<void*> (&time), 1, MPI_DOUBLE, 0, tmsh.comm);
    MPI_Bcast (static_cast<void*> (&savecount), 1, MPI_DOUBLE, 0, tmsh.comm);
    MPI_Barrier (tmsh.comm); // tmsh.comm = MPI_COMM_WORLD
    
    // Print current time
    if(rank==0)
    {
      std::cout << "TIME = " << time << ", dt = " << deltat << std::endl;
      full_time_vector.push_back (time);
    }
    
    
    // first step!
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.first_step(quadrant, deltat);
    }
    


    MPI_Allreduce(stp.sol_onehalf.data(), stp.local_extrema.data(), stp.sol_onehalf.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);

//    if (rank==0)
//    {
//      for (int ii=0;ii<stp.local_extrema.size();ii++)
//      std::cout << stp.local_extrema[ii] << " " << stp.sol_onehalf[ii] << std::endl;
//    }
//
//    if (rank==1)
//    {
//      for (int ii=0;ii<stp.local_extrema.size();ii++)
//      std::cout << stp.local_extrema[ii] << " " << stp.sol_onehalf[ii] << std::endl;
//    }
//    MPI_Barrier (MPI_COMM_WORLD);
//    MPI_Finalize ();
//    return 0;
    
    
//    for (int ii = 0; ii < gn_elements; ++ii)
//    {
//      if (stp.sol_onehalf[ii] < 0.)
//      {
//        stp.sol_onehalf[ii] = 0.;
//      }
//    }
    
    
    
    // second step!
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
    {
      stp.second_step(quadrant, incr);
    }

    
    incr.assemble ();
    TOC("Compute step");
    

    
    
    
    TIC();
    for (auto kk = 0; kk < incr.get_owned_data ().size (); kk++)
    {
      // this is the increment, sol.get_owned_data () is a vector probably because it has .assign function
      sol.get_owned_data ()[kk] += deltat * incr.get_owned_data ()[kk] / mass.get_owned_data ()[kk];
    }
    sol.assemble (replace_op);
    
//    for (auto quadrant = tmsh.begin_quadrant_sweep ();
//         quadrant != tmsh.end_quadrant_sweep (); ++quadrant)
//    {
//      for (int ii = 0; ii < 4; ++ii)
//      {
//        if (! quadrant->is_hanging (ii) )
//        {
//          auto & hdof = sol [ordh (quadrant->gt (ii) )];
//          hdof = hdof>0. ? hdof : 0.;
//        }
//      }
//    }
    
    TOC("Apply increment");
    
    // Save solution
    if (savecount >= SAVEDT) {
      TIC();
      if (rank == 0)
        std::cout << "savecount = " << savecount << std::endl;
      count++;
      save_time_vector.push_back (time);
      sprintf(filename, "results/swe_h_%4.4d",   count);
      tmsh.octbin_export (filename, sol, ordh);
      sprintf(filename, "results/swe_Ux_%4.4d",  count);
      tmsh.octbin_export (filename, sol, ordUx);
      sprintf(filename, "results/swe_Uy_%4.4d",  count);
      tmsh.octbin_export (filename, sol, ordUy);
      sprintf(filename, "results/swe_Z_%4.4d",  count);
      tmsh.octbin_export (filename, Z);
      savecount = 0.0;
      TOC("Exporting solution");
      
    }
    
    
  }
  
  
  if (rank == 0)
  {
    ColumnVector vtmp (save_time_vector.size ());
    std::copy (save_time_vector.begin (), save_time_vector.end (), vtmp.fortran_vec ());
    octave_io_mode m;
    octave_io_open ("timesteps.octbin", gz_write_mode, &m);
    octave_save ("save_time", vtmp);
    vtmp.resize (full_time_vector.size ());
    std::copy (full_time_vector.begin (), full_time_vector.end (), vtmp.fortran_vec ());
    octave_save ("full_time", vtmp);
    octave_io_close ();
  }
  
  
  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }
  
  MPI_Finalize ();
  
  
  return 0;
  
}
