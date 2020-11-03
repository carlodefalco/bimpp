#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators_3d.h>
#include <tmesh_3d.h>
#include <octave_file_io.h>

#include <vector>
#include <cassert>

int
write_example_connectivity3 (const char* filename)
{
  std::vector<double> p =
    {0.00, 0.00, 0.00, 0.00,
     0.25, 0.25, 0.25, 0.25,
     0.75, 0.75, 0.75, 0.75,
     1.00, 1.00, 1.00, 1.00,
     0.00, 0.00, 0.00, 0.00,
     0.25, 0.25, 0.25, 0.25,
     0.75, 0.75, 0.75, 0.75,
     1.00, 1.00, 1.00, 1.00,
     0.00, 0.00, 0.00, 0.00,
     0.25, 0.25, 0.25, 0.25,
     0.75, 0.75, 0.75, 0.75,
     1.00, 1.00, 1.00, 1.00,
     0.00, 0.00, 0.00, 0.00,
     0.25, 0.25, 0.25, 0.25,
     0.75, 0.75, 0.75, 0.75,
     1.00, 1.00, 1.00, 1.00,
     //
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     0.00, 0.25, 0.75, 1.00,
     //
     0.00, 0.00, 0.00, 0.00,
     0.00, 0.00, 0.00, 0.00,
     0.00, 0.00, 0.00, 0.00,
     0.00, 0.00, 0.00, 0.00,
     0.25, 0.25, 0.25, 0.25,
     0.25, 0.25, 0.25, 0.25,
     0.25, 0.25, 0.25, 0.25,
     0.25, 0.25, 0.25, 0.25,
     0.75, 0.75, 0.75, 0.75,
     0.75, 0.75, 0.75, 0.75,
     0.75, 0.75, 0.75, 0.75,
     0.75, 0.75, 0.75, 0.75,
     1.00, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00};
  
  std::vector<octave_idx_type> t =
    { 1,  5,  2,  6, 17, 21, 18, 22,  1,
      2,  6,  3,  7, 18, 22, 19, 23,  1,
      3,  7,  4,  8, 19, 23, 20, 24,  1,
      5,  9,  6, 10, 21, 25, 22, 26,  1,
      6, 10,  7, 11, 22, 26, 23, 27,  1,
      7, 11,  8, 12, 23, 27, 24, 28,  1,
      9, 13, 10, 14, 25, 29, 26, 30,  1,
      10, 14, 11, 15, 26, 30, 27, 31,  1,
      11, 15, 12, 16, 27, 31, 28, 32,  1,
      17, 21, 18, 22, 33, 37, 34, 38,  1,
      18, 22, 19, 23, 34, 38, 35, 39,  1,
      19, 23, 20, 24, 35, 39, 36, 40,  1,
      21, 25, 22, 26, 37, 41, 38, 42,  1,
      22, 26, 23, 27, 38, 42, 39, 43,  1,
      23, 27, 24, 28, 39, 43, 40, 44,  1,
      25, 29, 26, 30, 41, 45, 42, 46,  1,
      26, 30, 27, 31, 42, 46, 43, 47,  1,
      27, 31, 28, 32, 43, 47, 44, 48,  1,
      33, 37, 34, 38, 49, 53, 50, 54,  1,
      34, 38, 35, 39, 50, 54, 51, 55,  1,
      35, 39, 36, 40, 51, 55, 52, 56,  1,
      37, 41, 38, 42, 53, 57, 54, 58,  1,
      38, 42, 39, 43, 54, 58, 55, 59,  1,
      39, 43, 40, 44, 55, 59, 56, 60,  1,
      41, 45, 42, 46, 57, 61, 58, 62,  1,
      42, 46, 43, 47, 58, 62, 59, 63,  1,
      43, 47, 44, 48, 59, 63, 60, 64,  1};
  
  // save data to file
  Matrix oct_p (p.size() / 3, 3, 0.0);
  Array<octave_idx_type> oct_t (dim_vector (9, t.size() / 9), 0);
  
  std::copy_n (p.begin (), p.size (), oct_p.fortran_vec ());
  oct_p = oct_p.transpose ();
  std::copy_n (t.begin (), t.size (), oct_t.fortran_vec ());
  
  octave_scalar_map the_map;
  the_map.assign ("p", oct_p);
  the_map.assign ("t", oct_t);
  
  octave_io_mode m = gz_write_mode;
  int CHK;
  CHK = octave_io_open (filename, m, &m); assert (CHK == 0);
  CHK = octave_save ("msh", octave_value (the_map)); assert (CHK == 0);
  CHK = octave_io_close (); assert (CHK == 0);
  
  return 0;
}

static int
uniform_refinement (tmesh_3d::quadrant_iterator quadrant)
{ return 1; }

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;
  
  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity3 ("p4est_operator_test_3d.octbin");

  tmsh.read_connectivity ("p4est_operator_test_3d.octbin.gz");
  
  // Uniform refinement.
  recursive = 0; partforcoarsen = 1;
  for (int cycle = 0; cycle < 4; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  tmsh.vtk_export ("p4est_operator_test_3d");
  
  // Assemble advection-diffusion matrix.
  sparse_matrix A;
  A.resize(tmsh.num_global_nodes());
  
  std::vector<double> alpha(tmsh.num_local_quadrants (), 1);
  std::vector<double> psi(tmsh.num_global_nodes (), 0);
  
  double x, y, z, rho;
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 8; ++ii)
        {
          if (! quadrant->is_hanging (ii))
            {
              x = quadrant->p(0, ii);
              y = quadrant->p(1, ii);
              z = quadrant->p(2, ii);
              
              rho = std::sqrt(x * x + y * y + z * z);
              
              if (rho >= 0.8 && rho <= 0.9)
                psi[quadrant->gt(ii)] = -(2 * rho - 0.8) / 1e-2;
              else if (rho >= 0.9)
                psi[quadrant->gt(ii)] = -0.2 / 1e-2;
            }
        }
    }
  
  // Reduce coefficients.
  std::vector<double> global_psi(tmsh.num_global_nodes(), 0);
  MPI_Allreduce(psi.data(), global_psi.data(), psi.size(),
                MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  bim3a_advection_diffusion (tmsh, alpha, global_psi, A);
  
  // Assemble right-hand side.
  std::vector<double> rhs(tmsh.num_global_nodes (), 0);
  
  std::vector<double> f(tmsh.num_local_quadrants (), 0);
  std::vector<double> g(tmsh.num_global_nodes (), 0);
  
  bim3a_rhs (tmsh, f, g, rhs);
  
  // Set boundary conditions.
  dirichlet_bcs3 bcs;
  bcs.push_back (std::make_tuple(0,  0,
                  [] (double x, double y, double z) { return 0.3; }));
  bcs.push_back (std::make_tuple(0,  2,
                  [] (double x, double y, double z) { return 0.3; }));
  bcs.push_back (std::make_tuple(0,  4,
                  [] (double x, double y, double z) { return 0.3; }));
  bcs.push_back (std::make_tuple(26, 1,
                  [] (double x, double y, double z) { return 0; }));
  bcs.push_back (std::make_tuple(26, 3,
                  [] (double x, double y, double z) { return 0; }));
  bcs.push_back (std::make_tuple(26, 5,
                  [] (double x, double y, double z) { return 0; }));
  
  bim3a_dirichlet_bc (tmsh, bcs, A, rhs);
  
  // Solve problem.
  std::cout << "Solving linear system." << std::endl;
  
  mumps mumps_solver;
  
  std::vector<double> vals;
  std::vector<int> irow, jcol;
  
  A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
  
  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);
  
  // Reduce rhs (so that rank 0 has the actual rhs).
  std::vector<double> global_rhs(tmsh.num_global_nodes(), 0);
  MPI_Allreduce(rhs.data(), global_rhs.data(), rhs.size(),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (rank == 0)
    mumps_solver.set_rhs (global_rhs);
  
  // Solve.
  mumps_solver.analyze ();
  mumps_solver.factorize ();
  mumps_solver.solve ();
  mumps_solver.cleanup ();
  
  // Export solution.
  MPI_Bcast(global_rhs.data(), global_rhs.size(),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
  tmsh.octbin_export ("p4est_operator_test_3d_output", global_rhs);
  
  MPI_Finalize ();
  
  return 0;
}
