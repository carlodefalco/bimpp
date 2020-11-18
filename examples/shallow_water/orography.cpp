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

static constexpr char LOADFILENAME_1[255] = "orography.octbin.gz";
static constexpr char LOADFILENAME_2[255] = "basin.octbin.gz";
static constexpr char SAVEFILENAME_1[255] = "orography_tmsh";
static constexpr char SAVEFILENAME_2[255] = "orography_tmsh";
static constexpr char VARNAME_1[255] = "dem";
static constexpr char VARNAME_2[255] = "basin";

static constexpr double res = 5;
static constexpr double Nx = 1998;
static constexpr double Ny = 1829;
static constexpr double L = res*(Nx-1);
static constexpr double H = res*(Ny-1);
static std::vector<double>   dem;
static std::vector<double>   basin_mask;
static constexpr int NUM_REFINEMENTS = 9; // 3
static constexpr int NUM_TREFINEMENTS = 3;

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

//static int
//refine_function (tmesh::quadrant_iterator quadrant)  // -1, 0, 1
//{
//  static double x[4];
//  static double y[4];
//  static double xm;
//  static double ym;
//
//  // x, y coord. of cell center
//  xm = quadrant->centroid (0);
//  ym = quadrant->centroid (1);
//
//  for (int ii = 0; ii < 4; ++ii) {
//    x[ii] = quadrant->p (0, ii);
//    y[ii] = quadrant->p (0, ii);
//  }
//
//
//  return (quadrant->centroid (0) > 500 ? 1 : 0);
//
//}

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
  
  ii = i_y + Ny*i_x;
  
  return(std::array<int,3>{{ ii,int(i_x),int(i_y) }});
}

static int
raster_2_vector(const double& i_x,
                const double& i_y)
{
  static int ii;
  ii = i_y + Ny*i_x;
  
  return(ii);
}



static int
refine_function (tmesh::quadrant_iterator quadrant)  // -1, 0, 1
{
  
  static double xm;
  static double ym;
  
  static std::array<int,3> ids;
  
  static int i_x_plus;
  static int i_x_minus;
  static int i_y_plus;
  static int i_y_minus;
  
  // x, y coord. of cell center
  xm = quadrant->centroid (0);
  ym = quadrant->centroid (1);
  
  ids = global_coord_2_raster(xm,ym);
  
  const auto & id_m = ids[0];
  
  const auto & i_x_m = ids[1];
  const auto & i_y_m = ids[2];
  
  
  i_x_plus = i_x_m+1;
  i_x_plus = i_x_plus >= Nx ? (Nx-1) : i_x_plus;
  
  i_x_minus = i_x_m-1;
  i_x_minus = i_x_minus < 0 ? 0 : i_x_minus;
  
  i_y_plus = i_y_m+1;
  i_y_plus = i_y_plus >= Ny ? (Ny-1) : i_y_plus;
  
  i_y_minus = i_y_m-1;
  i_y_minus = i_y_minus < 0 ? 0 : i_y_minus;
  
  return(
  basin_mask[id_m] ?
  !(basin_mask[raster_2_vector(i_x_plus,i_y_m)] &&
  basin_mask[raster_2_vector(i_x_minus,i_y_m)] &&
  basin_mask[raster_2_vector(i_x_m,i_y_plus)] &&
  basin_mask[raster_2_vector(i_x_m,i_y_minus)] &&
  basin_mask[raster_2_vector(i_x_minus,i_y_minus)] &&
  basin_mask[raster_2_vector(i_x_plus,i_y_plus)] &&
  basin_mask[raster_2_vector(i_x_minus,i_y_plus)] &&
  basin_mask[raster_2_vector(i_x_plus,i_y_minus)]) ? 1 : 0 : 0);
  
//  return(
//         basin_mask[id_m] ?
//         !(basin_mask[raster_2_vector(i_x_plus,i_y_m)] &&
//           basin_mask[raster_2_vector(i_x_minus,i_y_m)] &&
//           basin_mask[raster_2_vector(i_x_m,i_y_plus)] &&
//           basin_mask[raster_2_vector(i_x_m,i_y_minus)]) ? 1 : 0 : 0);
  
}






// Re-Define tic and toc to add an MPI_Barrier
#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }

using Q1  = q1_vec<distributed_vector>;  // Typedef for distributed q_1 vector


int
main (int argc, char **argv)
{
  // Manegement of solutions ordering
  ordering // bim_ordering.h
  ord0 = [] (tmesh::idx_t gt) -> size_t { return dof_ordering<1, 0> (gt); };
  
  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  
  // Buffer for export filename
  char filename[255]="";
  
  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
  
  tmesh::idx_t ln_nodes = tmsh.num_owned_nodes ();
  
    /// Allocate initial data container
  Q1 sol (ln_nodes);
  sol.get_owned_data ().assign (sol.get_owned_data ().size (), 0.0);

  TIC();
  octave_io_mode m_in = gz_read_mode, m_out = gz_read_mode;
  octave_value v;
  
  octave_io_open (LOADFILENAME_1, m_in, &m_out);
  octave_load (VARNAME_1, v);
  Matrix M = v.matrix_value ();
  dem.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), dem.begin ());
  
  octave_io_open (LOADFILENAME_2, m_in, &m_out);
  octave_load (VARNAME_2, v);
  M = v.matrix_value ();
  basin_mask.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), basin_mask.begin ());
  TOC("Load data matrix");
  
  
  
  TIC ();
  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);
  TOC ("Uniform refinement");
  
  // Initialize initial datas
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
  {
    
    for (int ii = 0; ii < 4; ++ii)
    {
      if (! quadrant->is_hanging (ii)){
        double xx=quadrant->p(0,ii);
        double yy=quadrant->p(1,ii);
//        std::cout << get_raster_ref(xx, yy) << std::endl;
        sol[ord0(quadrant->gt (ii))] = dem[global_coord_2_raster(xx,yy)[0]];
      }
      
      else
      {
        sol[ord0(quadrant->gparent(0,ii))] +=0.;
        sol[ord0(quadrant->gparent(1,ii))] +=0.;
      }
    }
  }

  TIC ();
  /// Save initial conditions
  sprintf(filename, "initial_orography");
  tmsh.octbin_export (filename, sol, ord0);
  TOC("Save initial orography");
  
  
  TIC ();
  for (int ii = 0; ii < NUM_TREFINEMENTS; ++ii) {
    tmsh.set_refine_marker (refine_function);
    tmsh.refine ();
  }
  TOC ("Non-uniform refinement");
  
  TIC ();
  /// Save the p4est and connectivity to a file.
  tmsh.save (SAVEFILENAME_1);
  tmsh.vtk_export (SAVEFILENAME_1);
  TOC ("Save");
  return 0;


  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  MPI_Finalize ();

  return 0;
}
