#include <bim_timing.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
#include <cassert>

char filename[255];

static int
bottom_refinement (tmesh::quadrant_iterator quadrant)
{
  double ycoord;
  double top = std::numeric_limits<double>::lowest ();
  for (int ii = 0; ii < 4; ++ii)
    {
      ycoord = quadrant->p(1, ii);
      top  = top < ycoord ? ycoord : top;
    }
  return ((top <= 0.5) ? 1 : 0);
}

static int
left_refinement (tmesh::quadrant_iterator quadrant)
{
  double xcoord;
  double right = std::numeric_limits<double>::lowest ();
  for (int ii = 0; ii < 4; ++ii)
    {
      xcoord = quadrant->p(0, ii);
      right  = right < xcoord ? xcoord : right;
    }
  return ((right <= 0.5) ? 1 : 0);
}

static int
uniform_refinement (tmesh::quadrant_iterator quadrant)
{ return 1; }

int
main (int argc, char **argv)
{

  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh                 tmsh;
  double                xcoord = 0.0;
  double                ycoord = 0.0;

  MPI_Init (&argc, &argv);

  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_conn_test_hanging.octbin.gz");

  tmsh.read_connectivity ("p4est_conn_test_hanging.octbin.gz");

  tmsh.set_refine_marker (uniform_refinement);
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.set_refine_marker (bottom_refinement);
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.set_refine_marker (left_refinement);
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_conn_test_hanging");

  std::vector<double> p(2 * tmsh.num_local_nodes());
  Cell oct_t(4, tmsh.num_local_elems());
  ColumnVector parents(2, 0);
  
  std::array<int, 4> local_idx = {0, 1, 3, 2};
  
  int row = 0;
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      row = 0;
      for (auto ii : local_idx)
        {
          if (! quadrant->is_hanging(ii))
            {
              xcoord = quadrant->p(0, ii);
              ycoord = quadrant->p(1, ii);
              
              oct_t(row++, quadrant->get_forest_quad_idx()) = octave_value(quadrant->t(ii));
              
              p[2 * quadrant->t(ii) + 0] = xcoord;
              p[2 * quadrant->t(ii) + 1] = ycoord;
            }
          else
            {
              parents(0) = quadrant->parent(0, ii);
              parents(1) = quadrant->parent(1, ii);
              
              oct_t(row++, quadrant->get_forest_quad_idx()) = octave_value(parents);
            }
        }
    }
  
  Matrix oct_p(2, p.size() / 2, 0.0);
  Array<int> oct_children (dim_vector(4, tmsh.num_local_elems()), 0);
  
  std::copy_n (p.begin (), p.size (), oct_p.fortran_vec ());
  
  octave_scalar_map the_map;
  the_map.assign ("p", oct_p);
  the_map.assign ("t", oct_t);
  the_map.assign ("children", oct_children);
  
  octave_io_mode m = gz_write_mode;
  sprintf(filename, "p4est_conn_test_hanging_output_%4.4d.octbin.gz", rank);
  assert (octave_io_open (filename, m, &m) == 0);
  assert (octave_save ("msh", octave_value (the_map)) == 0);
  assert (octave_io_close () == 0);
      
  MPI_Finalize ();
  return 0;
}
