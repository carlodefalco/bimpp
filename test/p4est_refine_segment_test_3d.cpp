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
  Matrix oct_p (8, 3, 0);
  Array<int> oct_t (dim_vector(9, 1), 0);
    
  std::copy_n (p.begin (), p.size (), oct_p.fortran_vec ());
  oct_p = oct_p.transpose ();
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
    
  for (int segment = 0; segment < segment_list.size(); ++segment)
    {
      Point A = segment_list[segment][0];
      Point B = segment_list[segment][1];
        
      if (point_in_cube (A, lbf, rtb) || point_in_cube (B, lbf, rtb))
        return 1;
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
  s_lst.push_back ({Point({0.1, 0.1, 0.1}), Point({0.5, 0.5, 0.1})});

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
  for (int cycle = 0; cycle < 16; ++cycle)
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
