#include <cassert>
#include <vector>

#include <bim_timing.h>
#include <octave_file_io.h>
#include <tmesh.h>

#include <mosfet_connectivity_2d.h>
#include <mosfet_doping_2d.h>



char filename[255];

static int
doping_driven_refinement (tmesh::quadrant_iterator quadrant)
{

  constexpr double L = 3.0e-6;
  constexpr double H = 1.0e-5;

  double maxy = 0, miny = 0, y = 0;
  double xcoord, ycoord;
  
  double top = std::numeric_limits<double>::lowest ();

  maxy = std::numeric_limits<double>::lowest ();
  miny = std::numeric_limits<double>::max ();

  for (int ii = 0; ii < 4; ++ii)
    {

      xcoord = quadrant->p(0, ii);
      ycoord = quadrant->p(1, ii);
      
      y = signedlog (doping (xcoord, ycoord, L, H));

      maxy = maxy < y ? y : maxy;
      miny = miny > y ? y : miny;
      top  = top < ycoord ? ycoord : top;
    }

  double delta = maxy - miny;
  return ((top <= 0 && delta > .1) ? 1 : 0);

}

static int
coarsen_right_half (tmesh::quadrant_iterator quadrant)
{
  constexpr double L = 3.0e-6;
  constexpr double H = 1.0e-5;

  double minx = std::numeric_limits<double>::max ();
  double xcoord, ycoord;
  
  for (tmesh::idx_t jj = 0; jj < 4; ++jj)
    {
      xcoord = quadrant->p(0, jj);
      minx = minx > xcoord ? xcoord : minx;
    }
  
  return (minx > L/2.0 ? 1 : 0);
  
}

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_conn_test.octbin.gz");

  tmsh.read_connectivity ("p4est_conn_test.octbin.gz");
  
  recursive = 0;
  partforcoarsen = 0;
  
  tmsh.set_refine_marker (doping_driven_refinement);  
  tmsh.refine (recursive, partforcoarsen);
  tmsh.set_refine_marker (doping_driven_refinement);
  tmsh.refine (recursive, partforcoarsen);

  /*recursive = 1;
  tmsh.set_coarsen_marker (coarsen_right_half);*
  tmsh.coarsen (recursive, partforcoarsen);*/
  
  tmsh.vtk_export ("p4est_conn_test");
  
  double xcoord = 0.0;
  double ycoord = 0.0;
  
  std::vector<double> p(2 * tmsh.num_local_nodes());
  std::vector<double> p_hanging;
  std::vector<octave_idx_type> t;
  
  std::array<int, 4> local_idx = {0, 1, 3, 2};
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (auto ii : local_idx)
        {
          if (! quadrant->is_hanging(ii))
            {
              xcoord = quadrant->p(0, ii);
              ycoord = quadrant->p(1, ii);
              
              t.push_back(quadrant->t(ii));
              
              p[2 * t.back() + 0] = xcoord;
              p[2 * t.back() + 1] = ycoord;
            }
          else
            {
              switch (ii)
                {
                  case 0:
                    xcoord = quadrant->p(0, 2);
                    ycoord = quadrant->p(1, 1);
                    break;
                    
                  case 1:
                    xcoord = quadrant->p(0, 3);
                    ycoord = quadrant->p(1, 0);
                    break;
                    
                  case 2:
                    xcoord = quadrant->p(0, 0);
                    ycoord = quadrant->p(1, 3);
                    break;
                    
                  case 3:
                    xcoord = quadrant->p(0, 1);
                    ycoord = quadrant->p(1, 2);
                    break;
                }
                
                p_hanging.push_back(xcoord);
                p_hanging.push_back(ycoord);
                
                t.push_back(p.size() / 2 + p_hanging.size() / 2 - 1);
            }
        }
    }
  
  Matrix oct_p(2, p.size() / 2, 0.0);
  Matrix oct_p_hanging(2, p_hanging.size() / 2, 0.0);
  
  Array<octave_idx_type> oct_t(dim_vector(4, t.size() / 4), 0);
  Array<octave_idx_type> oct_children (oct_t);
  
  std::copy_n (p.begin (), p.size (), oct_p.fortran_vec ());
  std::copy_n (p_hanging.begin (), p_hanging.size (), oct_p_hanging.fortran_vec ());
  std::copy_n (t.begin (), t.size (), oct_t.fortran_vec ());
  
  octave_scalar_map the_map;
  the_map.assign ("p", oct_p);
  the_map.assign ("p_hanging", oct_p_hanging);
  the_map.assign ("t", oct_t);
  the_map.assign ("children", oct_children);
  
  octave_io_mode m = gz_write_mode;
  sprintf(filename, "p4est_conn_test_output_%4.4d.octbin.gz", rank);
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_save ("msh", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);
  
  MPI_Finalize ();
  return 0;
}
