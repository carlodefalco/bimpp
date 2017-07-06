#include <simple_connectivity_2d.h>
#include <octave_file_io.h>

int
write_example_connectivity (const char* filename)
{
  std::vector<double> p = {0., 1., 1., 0.,
                           0., 0., 1., 1.};

  std::vector<int> t = {1, 2, 3, 4, 1};

  // save data to file
  Matrix oct_p (4, 2, 0.0);
  Array<int> oct_t (dim_vector (5, 1), 0);

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
