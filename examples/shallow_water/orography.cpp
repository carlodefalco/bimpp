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

static constexpr char LOADFILENAME[255] = "orography_tiff";
static constexpr char SAVEFILENAME[255] = "orography_tmsh";
static constexpr char VARNAME[255] = "orography";

static constexpr double L = 1000;
static constexpr double H = 1000;
static std::vector<double>   data_matrix;
static constexpr int NUM_REFINEMENTS = 3;
static constexpr int NUM_TREFINEMENTS = 9;

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
refine_function (tmesh::quadrant_iterator quadrant)
{
  static double x[4];
  static double y[4];
  static double xm;
  static double ym;

  xm = quadrant->centroid (0);
  ym = quadrant->centroid (1);
  
  for (int ii = 0; ii < 4; ++ii) {
    x[ii] = quadrant->p (0, ii);
    y[ii] = quadrant->p (0, ii);
  }

  return (quadrant->centroid (0) > 0.0 ? 1 : 0);

}

// Re-Define tic and toc to add an MPI_Barrier
#define TIC() MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
#define TOC(S) MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { toc (S); }


int
main (int argc, char **argv)
{

  // Initialize MPI
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  TIC();
  octave_io_mode m_in = gz_read_mode, m_out = gz_read_mode;
  octave_io_open (LOADFILENAME, m_in, &m_out);
  octave_value v;
  octave_load (VARNAME, v);
  Matrix M = v.matrix_value ();
  data_matrix.resize (M.numel ());
  std::copy (M.fortran_vec (), M.fortran_vec () + M.numel (), data_matrix.begin ());
  TOC("Load data matrix");
  
  /// Generate the mesh in 2d
  tmesh tmsh;
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  TIC ();
  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);
  TOC ("Uniform refinement");

  TIC ();
  for (int ii = 0; ii < NUM_TREFINEMENTS; ++ii) {
    tmsh.set_refine_marker (refine_function);
    tmsh.refine ();
  }
  TOC ("Non-uniform refinement");

  TIC ();
  /// Save the p4est and connectivity to a file.
  tmsh.save (SAVEFILENAME);
  TOC ("Save");


  // Close MPI and print report
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

  MPI_Finalize ();

  return 0;
}
