#include <tmesh.h>
#include <bim_timing.h>
#include <octave_file_io.h>
#include <vector>
#include <cstdio>
char name[255];
char step[255];

using Point = std::array<double, 2>;
using Segment = std::array<Point, 2>;

int
unit_square (const char * filename)
{
  std::vector<double> p = {0, 1, 1, 0,
                           0, 0, 1, 1};
    
  std::vector<int> t = {1, 2, 3, 4, 1};
    
  // Save data to file.
  Matrix oct_p (4, 2, 0);
#ifdef HAVE_OCTAVE_44
  Array<octave_int32> 
#else
    Array<int> 
#endif
    oct_t (dim_vector(5, 1), 0);
    
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
uniform_refinement (tmesh::quadrant_iterator quadrant)
{ return 1; }

static int
segment_list_refinement (tmesh::quadrant_iterator quadrant,
                         const std::vector<Segment> & segment_list)
{
  double x0 = quadrant->p (0, 0);
  double y0 = quadrant->p (1, 0);

  double x1 = quadrant->p (0, 1);
  double y1 = quadrant->p (1, 2);
    
  double t = 0;
    
  for (int segment = 0; segment < segment_list.size(); ++segment)
    {
      Point A = segment_list[segment][0];
      Point B = segment_list[segment][1];
        
      double xA = A[0]; double yA = A[1];
      double xB = B[0]; double yB = B[1];
        
      // Ignore segments that don't intersect current quadrant.
      if ((xA <= x0 && xB <= x0)
          || (xA >= x1 && xB >= x1)
          || (yA <= y0 && yB <= y0)
          || (yA >= y1 && yB >= y1))
        { continue; }
        
      /*
       * Sides numbering:
       * 
       *        3
       *   +---------+
       *   |         |
       * 0 |         | 1
       *   |         |
       *   +---------+
       *        2
       */
        
      // If AB is not horizontal.
      if (yA != yB)
        {
          // Side 2.
          t = (y0 - yB) / (yA - yB);
            
          double x_intersect = xA * t + xB * (1 - t);
            
          if (x_intersect >= x0 && x_intersect <= x1)
            return true;            
            
          // Side 3.
          t = (y1 - yB) / (yA - yB);
            
          x_intersect = xA * t + xB * (1 - t);
            
          if (x_intersect >= x0 && x_intersect <= x1)
            return true;
        }
        
      // If AB is not vertical.
      if (xA != xB)
        {
          // Side 0.
          t = (x0 - xB) / (xA - xB);
            
          double y_intersect = yA * t + yB * (1 - t);
            
          if (y_intersect >= y0 && y_intersect <= y1)
            return true;
            
          // Side 1.
          t = (x1 - xB) / (xA - xB);
            
          y_intersect = yA * t + yB * (1 - t);
            
          if (y_intersect >= y0 && y_intersect <= y1)
            return true;
        }
    }
    
  return false;
}

int main(int argc, char ** argv)
{
  MPI_Init (&argc, &argv);
  
  int      recursive, partforcoarsen, balance;
  MPI_Comm mpicomm = MPI_COMM_WORLD;  
  int      rank, size;
  tmesh    tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);
    
  MPI_Barrier (MPI_COMM_WORLD); { if (rank == 0) tic (); }
  
  // Create mesh.
  if (rank == 0)
    unit_square ("p4est_unitsquare.octbin.gz");    
  tmsh.read_connectivity ("p4est_unitsquare.octbin.gz");

  // Define marking for adaptive refinement.
  std::vector<Segment> segment_list;
  segment_list.push_back ({Point({0.10, 0.10}), Point({0.30, 0.30})});
  segment_list.push_back ({Point({0.30, 0.30}), Point({0.40, 0.90})});
  segment_list.push_back ({Point({0.25, 0.05}), Point({0.31, 0.29})});
  segment_list.push_back ({Point({0.31, 0.29}), Point({0.95, 0.45})});
  segment_list.push_back ({Point({0.57, 0.86}), Point({0.26, 0.25})});
  segment_list.push_back ({Point({0.88, 0.91}), Point({0.23, 0.38})});
  segment_list.push_back ({Point({0.07, 0.47}), Point({0.95, 0.36})});
  segment_list.push_back ({Point({0.71, 0.78}), Point({0.64, 0.46})});
  segment_list.push_back ({Point({0.82, 0.78}), Point({0.24, 0.57})});
  segment_list.push_back ({Point({0.10, 0.10}), Point({0.10, 0.90})});
  segment_list.push_back ({Point({0.10, 0.10}), Point({0.90, 0.10})});

  std::function<int (tmesh::quadrant_iterator)> segment_refinement =
    [segment_list] (tmesh::quadrant_iterator qi)
    { return segment_list_refinement(qi, segment_list); };


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
             sprintf (step, "*** Refinement and balancing %3.3d ***", refine_number++);
             toc (step); 
          }

        MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
        sprintf(name,"p4est_refine_segment_test_%4.4d",file_number++);
        tmsh.vtk_export (name);
        if (rank == 0) { toc ("*** Export ***"); }
    }
  
  // Refine according to segment_list.  
  for (int cycle = 0; cycle < 16; ++cycle)
    {
      // Adaptive refinement.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
      
      tmsh.set_refine_marker (segment_refinement);
      tmsh.refine (recursive, partforcoarsen);
      if (rank == 0) 
          { 
             sprintf (step, "*** Refinement and balancing %3.3d ***", refine_number++);
             toc (step); 
          }

      // Export mesh.
      MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }      
      sprintf(name,"p4est_refine_segment_test_%4.4d",file_number++);
      tmsh.vtk_export (name);
      if (rank == 0) { toc ("*** Export ***"); }  
    }
        
  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }
  
  MPI_Finalize ();
  return 0;
}
