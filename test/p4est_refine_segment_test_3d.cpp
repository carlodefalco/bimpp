#include <tmesh_3d.h>
#include <bim_timing.h>
#include <octave_file_io.h>
#include <vector>
#include <cstdio>
char name[255];
char step[255];

using Point = std::array<double, 3>;
using Segment = std::array<Point, 2>;

bool
point_in_cube (const Point pt, const Point lbf, const Point rtb)
{
  return (pt[0] > lbf[0] && pt[0] < rtb[0] &&
          pt[1] > lbf[1] && pt[1] < rtb[1] &&
          pt[2] > lbf[2] && pt[2] < rtb[2]);
}

int
unit_cube (const char * filename)
{
  std::vector<double> p = {0., 0., 0.,  1., 0., 0.,
                           0., 1., 0.,  1., 1., 0.,
                           0., 0., 1.,  1., 0., 1.,
                           0., 1., 1.,  1., 1., 1.};
    
  std::vector<int> t = {1, 2, 3, 4, 5, 6, 7, 8, 1};
    
  // Save data to file.
  Matrix oct_p (3, 8, 0);
  Array<int> oct_t (dim_vector(9, 1), 0);
    
  std::copy_n (p.begin (), p.size (), oct_p.fortran_vec ());
  std::copy_n (t.begin (), t.size (), oct_t.fortran_vec ());
    
  octave_scalar_map the_map;
  the_map.assign ("p", oct_p);
  the_map.assign ("t", oct_t);
    
  octave_io_mode m = gz_write_mode;
    
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_save ("msh", octave_value(the_map)) == 0);
  assert (octave_io_close () == 0);
    
  return 0;
}

static int
uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
{ return 1; }

static int
segment_list_refinement (tmesh_3d::quadrant_iterator quadrant,
                         const std::vector<Segment> & segment_list)
{
  Point lbf = {quadrant->p (0, 0), //left bottom front corner
               quadrant->p (1, 0),
               quadrant->p (2, 0)};
               
  Point rtb = {quadrant->p (0, 7), //right top back corner
               quadrant->p (1, 7),
               quadrant->p (2, 7)};

  bool pext;
  double t, inter_1, inter_2;
  std::array <int, 3> dim_i = {0, 1, 2};
    
  for (int segment = 0; segment < segment_list.size(); ++segment)
    {
      Point A = segment_list[segment][0];
      Point B = segment_list[segment][1];
        
      if (point_in_cube (A, lbf, rtb) || point_in_cube (B, lbf, rtb))
        return 1;
      
      // Ignore segments that don't intersect current quadrant.
      pext = false;
      for (auto ii : dim_i)
        pext = pext || (A[ii] < lbf[ii] && B[ii] < lbf[ii]) ||
                       (A[ii] > rtb[ii] && B[ii] > rtb[ii]);
      
      if (pext)
        continue;
        
      /*
       * Faces numbering:
       *
       *  +-------------+
       *  |\            |\
       *  | \     \5\   | \
       *  |  \  |3|     |  \
       *  |   +-------------+
       *  | 0 |         | 1 |
       *  +---|---------+   |
       *   \  |     |2|  \  |
       *    \ |   \4\     \ |
       *     \|            \|
       *      +-------------+
       *
       * Edges numbering:
       *        3
       *  +-------------+
       *  |\            |\
       *  | \6          | \7
       *  |10\      2   |11\
       *  |   +-------------+
       *  |   | 1       |   |
       *  +---|---------+   |
       *   \  |8         \  |9
       *    \4|           \5|
       *     \|     0      \|
       *      +-------------+
       */

      for (int ii = 0; ii < 3; ++ii)
        if (A[ii] != B[ii])
          {
            // Side 0, 2, 4
            t = (lbf[ii] - B[ii]) / (A[ii] - B[ii]);
              
            inter_1 = A[dim_i[(ii + 1) % 3]] * t +
                      B[dim_i[(ii + 1) % 3]] * (1 - t);
            inter_2 = A[dim_i[(ii + 2) % 3]] * t +
                      B[dim_i[(ii + 2) % 3]] * (1 - t);
              
            if (inter_1 >= lbf[dim_i[(ii + 1) % 3]] &&
                inter_1 <= rtb[dim_i[(ii + 1) % 3]] &&
                inter_2 >= lbf[dim_i[(ii + 2) % 3]] &&
                inter_2 <= rtb[dim_i[(ii + 2) % 3]])
              return 1;
              
            // Side 1, 3, 5
            t = (rtb[ii] - B[ii]) / (A[ii] - B[ii]);
              
            inter_1 = A[dim_i[(ii + 1) % 3]] * t +
                      B[dim_i[(ii + 1) % 3]] * (1 - t);
            inter_2 = A[dim_i[(ii + 2) % 3]] * t +
                      B[dim_i[(ii + 2) % 3]] * (1 - t);
              
            if (inter_1 >= lbf[dim_i[(ii + 1) % 3]] &&
                inter_1 <= rtb[dim_i[(ii + 1) % 3]] &&
                inter_2 >= lbf[dim_i[(ii + 2) % 3]] &&
                inter_2 <= rtb[dim_i[(ii + 2) % 3]])
              return 1;
          }
    }
    
  return 0;
}

int main(int argc, char ** argv)
{
  MPI_Init (&argc, &argv);
  
  int      recursive, partforcoarsen, balance;
  MPI_Comm mpicomm = MPI_COMM_WORLD;  
  int      rank, size;
  tmesh_3d tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);
    
  MPI_Barrier (MPI_COMM_WORLD); { if (rank == 0) tic (); }
  
  // Create mesh.
  if (rank == 0)
    unit_cube ("p4est_unitcube.octbin.gz");    
  tmsh.read_connectivity ("p4est_unitcube.octbin.gz");

  // Define marking for adaptive refinement.
  std::vector<Segment> s_lst;
  s_lst.push_back ({Point({.1, .7, .95}), Point({.1, .9, .95})});//P
  s_lst.push_back ({Point({.1, .9, .95}), Point({.2, .9, .95})});
  s_lst.push_back ({Point({.2, .9, .95}), Point({.2, .8, .95})});
  s_lst.push_back ({Point({.2, .8, .95}), Point({.1, .8, .95})});
  s_lst.push_back ({Point({.3, .7, .95}), Point({.3, .9, .95})});//A
  s_lst.push_back ({Point({.3, .9, .95}), Point({.4, .9, .95})});
  s_lst.push_back ({Point({.4, .9, .95}), Point({.4, .7, .95})});
  s_lst.push_back ({Point({.3, .8, .95}), Point({.4, .8, .95})});
  s_lst.push_back ({Point({.6, .9, .95}), Point({.5, .9, .95})});//C
  s_lst.push_back ({Point({.5, .9, .95}), Point({.5, .7, .95})});
  s_lst.push_back ({Point({.5, .7, .95}), Point({.6, .7, .95})});
  s_lst.push_back ({Point({.7, .7, .95}), Point({.8, .7, .95})});//S
  s_lst.push_back ({Point({.8, .7, .95}), Point({.8, .8, .95})});
  s_lst.push_back ({Point({.8, .8, .95}), Point({.7, .8, .95})});
  s_lst.push_back ({Point({.7, .8, .95}), Point({.7, .9, .95})});
  s_lst.push_back ({Point({.7, .9, .95}), Point({.8, .9, .95})});
  s_lst.push_back ({Point({.1, .6, .95}), Point({.2, .6, .95})});//2
  s_lst.push_back ({Point({.2, .6, .95}), Point({.2, .5, .95})});
  s_lst.push_back ({Point({.2, .5, .95}), Point({.1, .5, .95})});
  s_lst.push_back ({Point({.1, .5, .95}), Point({.1, .4, .95})});
  s_lst.push_back ({Point({.1, .4, .95}), Point({.2, .4, .95})});
  s_lst.push_back ({Point({.3, .4, .95}), Point({.3, .6, .95})});//0
  s_lst.push_back ({Point({.3, .6, .95}), Point({.4, .6, .95})});
  s_lst.push_back ({Point({.4, .6, .95}), Point({.4, .4, .95})});
  s_lst.push_back ({Point({.4, .4, .95}), Point({.3, .4, .95})});
  s_lst.push_back ({Point({.5, .5, .95}), Point({.6, .6, .95})});//1
  s_lst.push_back ({Point({.6, .6, .95}), Point({.6, .4, .95})});
  s_lst.push_back ({Point({.7, .5, .95}), Point({.7, .6, .95})});//8
  s_lst.push_back ({Point({.7, .6, .95}), Point({.8, .6, .95})});
  s_lst.push_back ({Point({.8, .6, .95}), Point({.8, .5, .95})});
  s_lst.push_back ({Point({.8, .5, .95}), Point({.7, .5, .95})});
  s_lst.push_back ({Point({.7, .5, .95}), Point({.7, .4, .95})});
  s_lst.push_back ({Point({.7, .4, .95}), Point({.8, .4, .95})});
  s_lst.push_back ({Point({.8, .4, .95}), Point({.8, .5, .95})});
  s_lst.push_back ({Point({.1, .2, 1.0}), Point({.8, .2, .95})});//--

  std::function<int (tmesh_3d::quadrant_iterator)> segment_refin =
    [s_lst] (tmesh_3d::quadrant_iterator qi)
    { return segment_list_refinement(qi, s_lst); };


  if (rank == 0) { toc ("*** Initialization ***"); }


  // Uniform refinement.
  recursive = 0;
  partforcoarsen = 1;
  
  int file_number = 0, refine_number = 0;
  for (int cycle = 0; cycle < 2; ++cycle)
    {
        MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
        tmsh.set_refine_marker (uniform_refinement);
        tmsh.refine(recursive, partforcoarsen);        
        if (rank == 0) 
          { 
             sprintf (step, "*** Refinement and balancing %3.3d ***",
                      refine_number++);
             toc (step); 
          }

        MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
        sprintf(name,"p4est_refine_segment_test_3d_%4.4d",
                file_number++);
        tmsh.vtk_export (name);
        if (rank == 0) { toc ("*** Export ***"); }
    }
  
  // Refine according to segment_list.  
  for (int cycle = 0; cycle < 12; ++cycle)
    {
      // Adaptive refinement.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      
      tmsh.set_refine_marker (segment_refin);
      tmsh.refine (recursive, partforcoarsen);
      
      if (rank == 0) 
          { 
             sprintf (step, "*** Refinement and balancing %3.3d ***",
                      refine_number++);
             toc (step); 
          }

      // Export mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }      
      sprintf(name,"p4est_refine_segment_test_3d_%4.4d",file_number++);
      tmsh.vtk_export (name);
      if (rank == 0) { toc ("*** Export ***"); }  
    }

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }
  
  MPI_Finalize ();
  return 0;
}
