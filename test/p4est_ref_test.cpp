#include <bim_timing.h>
#include <tmesh.h>

#include <p4est_bits.h>
#include <p4est_ghost.h>
#include <p4est_lnodes.h>
#include <p4est_vtk.h>
#include <p4est_connectivity.h>
#include <vector>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <octave_file_io.h>


static double
gaussian (const double x, const double c, const double s)
{ return (exp (- std::pow (x - c, 2) / (2 * std::pow (s, 2)))); }

static inline double
bool1 (const double x, const double L)
{ return (((x > L / 3.0) && (x < 2.0 * L / 3.0)) ? 1.0 : 0.0); }

static inline double
bool2 (const double x, const double L)
{ return (x <= L / 3.0 ? 1.0 : 0.0); }

static inline double
bool3 (const double x, const double L)
{ return (x >= 2.0 * L / 3.0 ? 1.0 : 0.0); }

static inline double
bool4 (const double x, const double L)
{ return (x < L / 4.0 ? 1.0 : 0.0); }

static inline double
bool5 (const double x, const double L)
{ return (x > 3.0 * L / 4.0 ? 1.0 : 0.0); }

static inline double
bool6 (const double x, const double L)
{ return (x >= L / 4.0 ? 1.0 : 0.0); }

static inline double
bool7 (const double x, const double L)
{ return (x <= 3.0 * L / 4.0 ? 1.0 : 0.0); }

static double
doping  (const double vxyz[3],
         const double L, const double H)
{
  constexpr double N_plus  = 1.0e25;
  constexpr double N_minus = 1.0e24;
  constexpr double P_plus  = 1.0e25;
  constexpr double P_minus = 1.0e24;

  double x = vxyz[0];
  double y = vxyz[1];

  double Na  = P_minus;
  Na += P_plus  * gaussian (y, -H / 2.0, H / 10.0) *
    (bool1 (x, L) +
     gaussian (x, L / 3.0, L / 10.0) * bool2 (x, L) +
     gaussian (x, 2.0 * L / 3.0, L / 10.0) * bool3 (x, L));

  double Nd  = 0.0;
  Nd += N_plus  * gaussian (y, 0, H / 4.0) *
    (bool4 (x, L) + bool5 (x, L) +
     gaussian (x, L / 4.0, L / 10.0) * bool6 (x, L) +
     gaussian (x, 3.0 * L / 4.0, L / 10.0) * bool7 (x, L));

  return (Nd - Na);
}

static inline double
signedlog (double x)
{ return (asinh (x / 2.0) / log (10.0)); }

static int
doping_driven_refinement (const double vxyz[12])
{

  constexpr double L = 3.0e-6;
  constexpr double H = 1.0e-5;

  double maxy = 0, miny = 0, y = 0;
  double top = vxyz[1];

  maxy = miny = y = signedlog (doping (&(vxyz[0]), L, H));

  for (int ii = 1; ii < 4; ++ii)
    {
      y = signedlog (doping (&(vxyz[0]) + 3 * ii, L, H));
      maxy = maxy < y ? y : maxy;
      miny = miny > y ? y : miny;
      top  = top < vxyz[1 + 3 * ii] ? vxyz[1 + 3 * ii] : top;
    }

  double delta = maxy - miny;
  return ((top <= 0 && delta > .1) ? 1 : 0);
}


static int
refine_fn (p4est_t * p4est, p4est_topidx_t tt,
           p4est_quadrant_t * quadrant)
{
  p4est_quadrant_t node;
  int i;
  double vxyz[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (i = 0; i < P4EST_CHILDREN; ++i)
    {
      p4est_quadrant_corner_node (quadrant, i, &node);
      p4est_qcoord_to_vertex (p4est->connectivity, tt, node.x, node.y,
                              &(vxyz[3 * i]));
    }

  return (doping_driven_refinement (vxyz));

}

static int
write_example_connectivity (const char* filename);


int
main (int argc, char **argv)
{

  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;
  
  MPI_Init (&argc, &argv);

  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_ref_test.octbin.gz");

  if (rank == 0)
    tmsh.read_connectivity ("p4est_ref_test.octbin.gz");

  tmsh.bcast_connectivity (0, mpicomm);
  tmsh.p4est = p4est_new (mpicomm, tmsh.conn, 0, NULL, NULL);

  recursive = 0;
  partforcoarsen = 0;

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  for (int k = 0; k < 13 ; ++k)
    {
      p4est_refine (tmsh.p4est, recursive, refine_fn, NULL);
      p4est_balance (tmsh.p4est, P4EST_CONNECT_FACE, NULL);
      p4est_partition (tmsh.p4est, partforcoarsen, NULL);
    }

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("refinement and balancing"); }


  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { tic (); }

  p4est_vtk_write_file (tmsh.p4est, NULL, "p4est_ref_test");

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0)
    { toc ("IO"); }

  p4est_destroy (tmsh.p4est);
  p4est_connectivity_destroy (tmsh.conn);

  if (rank == 0)
    { print_timing_report (); }

  MPI_Finalize ();
  return 0;

}


static int
write_example_connectivity (const char* filename)
{
  std::vector<double> p = {0.0e+00,    0.0e+00,    0.0e+00,    0.0e+00,    0.0e+00,
                           3.75e-07,   3.75e-07,   3.75e-07,   3.75e-07,   3.75e-07,
                           7.5e-07,    7.5e-07,    7.5e-07,    7.5e-07,    7.5e-07,
                           7.5e-07,    1.1250e-06, 1.1250e-06, 1.1250e-06, 1.1250e-06,
                           1.1250e-06, 1.1250e-06, 1.5e-06,    1.5e-06,    1.5e-06,
                           1.5e-06,    1.5e-06,    1.5e-06,    1.8750e-06, 1.8750e-06,
                           1.8750e-06, 1.8750e-06, 1.8750e-06, 1.8750e-06, 2.25e-06,
                           2.25e-06,   2.25e-06,   2.25e-06,   2.25e-06,   2.25e-06,
                           2.6250e-06, 2.6250e-06, 2.6250e-06, 2.6250e-06, 2.6250e-06,
                           3.0e-06,    3.0e-06,    3.0e-06,    3.0e-06,    3.0e-06,
                           //
                           -8.0e-06, -6.0e-06, -4.0e-06, -2.0e-06,  0.0e+00,
                           -8.0e-06, -6.0e-06, -4.0e-06, -2.0e-06,  0.0e+00,
                           -8.0e-06, -6.0e-06, -4.0e-06, -2.0e-06,  0.0e+00,
                           2.0e-06,  -8.0e-06, -6.0e-06, -4.0e-06, -2.0e-06,
                           0.0e+00,   2.0e-06, -8.0e-06, -6.0e-06, -4.0e-06,
                           -2.0e-06,  0.0e+00,  2.0e-06, -8.0e-06, -6.0e-06,
                           -4.0e-06, -2.0e-06,  0.0e+00,  2.0e-06, -8.0e-06,
                           -6.0e-06, -4.0e-06, -2.0e-06,  0.0e+00,  2.0e-06,
                           -8.0e-06, -6.0e-06, -4.0e-06, -2.0e-06,  0.0e+00,
                           -8.0e-06, -6.0e-06, -4.0e-06, -2.0e-06, 0.0e+00};

  std::vector<int> t = {15, 21,  22, 16, 2,
                        21, 27,  28, 22, 2,
                        27, 33,  34, 28, 2,
                        33, 39,  40, 34, 2,
                        1,  6,   7,  2,  3,
                        2,  7,   8,  3,  3,
                        3,  8,   9,  4,  3,
                        4,  9,  10,  5,  3,
                        6, 11,  12,  7,  3,
                        7, 12,  13,  8,  3,
                        8, 13,  14,  9,  3,
                        9, 14,  15, 10,  3,
                        11, 17,  18, 12, 3,
                        12, 18,  19, 13, 3,
                        13, 19,  20, 14, 3,
                        14, 20,  21, 15, 3,
                        17, 23,  24, 18, 3,
                        18, 24,  25, 19, 3,
                        19, 25,  26, 20, 3,
                        20, 26,  27, 21, 3,
                        23, 29,  30, 24, 3,
                        24, 30,  31, 25, 3,
                        25, 31,  32, 26, 3,
                        26, 32,  33, 27, 3,
                        29, 35,  36, 30, 3,
                        30, 36,  37, 31, 3,
                        31, 37,  38, 32, 3,
                        32, 38,  39, 33, 3,
                        35, 41,  42, 36, 3,
                        36, 42,  43, 37, 3,
                        37, 43,  44, 38, 3,
                        38, 44,  45, 39, 3,
                        41, 46,  47, 42, 3,
                        42, 47,  48, 43, 3,
                        43, 48,  49, 44, 3,
                        44, 49,  50, 45, 3};

  // save data to file
  Matrix oct_p (50, 2, 0.0);
  Array<int> oct_t (dim_vector (5, 36), 0);

  std::copy_n (p.begin (), p.size (), oct_p.fortran_vec ());
  oct_p = oct_p.transpose ();
  std::copy_n (t.begin (), t.size (), oct_t.fortran_vec ());

  octave_scalar_map the_map;
  the_map.assign ("p", oct_p);
  the_map.assign ("t", oct_t);

  octave_io_mode m = gz_write_mode;
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_save ("msh", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);


  return 0;
}
