#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <vector>
#include <cassert>

int
write_example_connectivity (const char* filename)
{
  std::vector<double> p = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
                           0.25000, 0.25000, 0.25000, 0.25000, 0.25000,
                           0.50000, 0.50000, 0.50000, 0.50000, 0.50000,
                           0.75000, 0.75000, 0.75000, 0.75000, 0.75000,
                           1.00000, 1.00000, 1.00000, 1.00000, 1.00000,
                           //
                           0.00000, 0.25000, 0.50000, 0.75000, 1.00000,
                           0.00000, 0.25000, 0.50000, 0.75000, 1.00000,
                           0.00000, 0.25000, 0.50000, 0.75000, 1.00000,
                           0.00000, 0.25000, 0.50000, 0.75000, 1.00000,
                           0.00000, 0.25000, 0.50000, 0.75000, 1.00000};
  
  std::vector<int> t = {1,  6,  7,  2, 1,
                        2,  7,  8,  3, 1,
                        3,  8,  9,  4, 1,
                        4,  9, 10,  5, 1,
                        6, 11, 12,  7, 1,
                        7, 12, 13,  8, 1,
                        8, 13, 14,  9, 1,
                        9, 14, 15, 10, 1,
                        11, 16, 17, 12, 1,
                        12, 17, 18, 13, 1,
                        13, 18, 19, 14, 1,
                        14, 19, 20, 15, 1,
                        16, 21, 22, 17, 1,
                        17, 22, 23, 18, 1,
                        18, 23, 24, 19, 1,
                        19, 24, 25, 20, 1};
  
  // save data to file
  Matrix oct_p (25, 2, 0.0);
  Array<int> oct_t (dim_vector (5, 16), 0);
  
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
  
  MPI_Init (&argc, &argv);

  mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);

  if (rank == 0)
    write_example_connectivity ("p4est_operator_test.octbin.gz");

  tmsh.read_connectivity ("p4est_operator_test.octbin.gz");
  
  // Uniform refinement.
  tmsh.set_refine_marker (uniform_refinement);
  
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_operator_test");
  
  // Assemble advection-diffusion matrix.
  sparse_matrix A;
  A.resize(tmsh.num_owned_nodes());
  
  std::vector<double> alpha(tmsh.num_local_elems (), 1);
  std::vector<double> psi(tmsh.num_local_nodes (), 0);
  
  double x = 0, y = 0, rho = 0;
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          x = quadrant->p(0, ii);
          y = quadrant->p(1, ii);
          
          rho = std::sqrt(x * x + y * y);
          
          if (rho >= 0.8 && rho <= 0.9)
            psi[quadrant->t(ii)] = -(2 * rho - 0.8) / 1e-2;
          else if (rho >= 0.9)
            psi[quadrant->t(ii)] = -0.2 / 1e-2;
        }
    }
  
  bim2a_advection_diffusion (tmsh, alpha, psi, A);
  
  // Assemble right-hand side.
  std::vector<double> rhs(tmsh.num_local_nodes (), 0);
  
  std::vector<double> f(tmsh.num_local_elems (), 0);
  std::vector<double> g(tmsh.num_local_nodes (), 0);
  
  bim2a_rhs (tmsh, f, g, rhs);
  
  // Set boundary conditions.
  dirichlet_bcs bcs;
  bcs.push_back (std::make_tuple( 0, 0, [] (double x, double y) { return 0.3; }));
  bcs.push_back (std::make_tuple( 0, 2, [] (double x, double y) { return 0.3; }));
  bcs.push_back (std::make_tuple(15, 1, [] (double x, double y) { return 0; }));
  bcs.push_back (std::make_tuple(15, 3, [] (double x, double y) { return 0; }));
  
  bim2a_dirichlet_bc (tmsh, bcs, A, rhs);
  
  // Solve problem.
  std::vector<double> vals;
  std::vector<int> irow, jcol;
  
  A.aij(vals, irow, jcol, 1);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  std::cout << "Solving linear system." << std::endl;
  
  mumps mumps_solver;
  
  if (rank == 0)
    mumps_solver.set_lhs_structure(A.rows(), irow, jcol);
  
  mumps_solver.analyze();
  
  if (rank == 0)
    mumps_solver.set_lhs_data(vals);
  
  mumps_solver.factorize();
  
  if (rank == 0)
    mumps_solver.set_rhs(rhs);
  
  mumps_solver.solve();
  
  // Export solution.
  tmsh.octbin_export ("p4est_operator_test_output", rhs);
  
  MPI_Finalize ();
  
  return 0;
}
