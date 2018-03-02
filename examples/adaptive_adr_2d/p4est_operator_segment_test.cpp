#include <tmesh.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <simple_connectivity_2d.h>
#include <bim_timing.h>
#include <vector>

using Point = std::array<double, 2>;
using Segment = std::array<Point, 2>;

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
       * 4 |         | 2
       *   |         |
       *   +---------+
       *        1
       */
        
      // If AB is not horizontal.
      if (yA != yB)
        {
          // Side 1.
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
          // Side 2.
          t = (x0 - xB) / (xA - xB);
            
          double y_intersect = yA * t + yB * (1 - t);
            
          if (y_intersect >= y0 && y_intersect <= y1)
            return true;
            
          // Side 4.
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
  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  // Define marking for adaptive refinement.
  std::vector<Segment> segment_list;
  segment_list.push_back ({Point({0.10, 0.10}), Point({0.30, 0.30})});
  segment_list.push_back ({Point({0.30, 0.30}), Point({0.40, 0.90})});
  segment_list.push_back ({Point({0.25, 0.05}), Point({0.31, 0.29})});
  //segment_list.push_back ({Point({0.31, 0.29}), Point({0.95, 0.45})});
  segment_list.push_back ({Point({0.57, 0.86}), Point({0.26, 0.25})});
  segment_list.push_back ({Point({0.88, 1.00}), Point({0.23, 0.00})});
  //segment_list.push_back ({Point({0.07, 0.47}), Point({0.95, 0.36})});
  segment_list.push_back ({Point({0.71, 0.78}), Point({0.64, 0.46})});
  segment_list.push_back ({Point({0.82, 0.78}), Point({0.24, 0.57})});
  //segment_list.push_back ({Point({0.10, 0.10}), Point({0.10, 0.90})});
  //segment_list.push_back ({Point({0.10, 0.10}), Point({0.90, 0.10})});

  std::function<int (tmesh::quadrant_iterator)> segment_refinement =
    [segment_list] (tmesh::quadrant_iterator qi)
    { return segment_list_refinement(qi, segment_list); };
  
  // Uniform refinement.
  recursive = 0;
  partforcoarsen = 1;
  
  for (int cycle = 0; cycle < 3; ++cycle)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine(recursive, partforcoarsen);
    }
  
  // Refine according to segment_list.  
  for (int cycle = 0; cycle < 16; ++cycle)
    {
      tmsh.set_refine_marker (segment_refinement);
      tmsh.refine (recursive, partforcoarsen);
    }
  
  if (rank == 0) { toc ("*** Mesh creation and refinement ***"); }
  
  // Assemble matrix.
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  
  sparse_matrix A;
  A.resize(tmsh.num_global_nodes());
  
  std::vector<double> alpha(tmsh.num_local_quadrants (), 1e-2);
  std::vector<double> psi(tmsh.num_local_nodes (), 0);
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      if (segment_refinement(quadrant))
        {
          alpha[quadrant->get_forest_quad_idx()] = 100;
        }
    }
    
  bim2a_advection_diffusion (tmsh, alpha, psi, A);
  
  if (rank == 0) { toc ("*** Matrix assembly ***"); }
  
  // Assemble right-hand side.
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  
  std::vector<double> rhs(tmsh.num_global_nodes (), 0);
  
  std::vector<double> f(tmsh.num_local_quadrants (), 0);
  std::vector<double> g(tmsh.num_local_nodes (), 0);
  
  bim2a_rhs (tmsh, f, g, rhs);
  
  if (rank == 0) { toc ("*** Right-hand side assembly ***"); }
  
  // Set boundary conditions.
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  
  dirichlet_bcs bcs;
  bcs.push_back (std::make_tuple(0, 2, [] (double x, double y) { return 0; }));
  bcs.push_back (std::make_tuple(0, 3, [] (double x, double y) { return 1; }));
  
  bim2a_dirichlet_bc (tmsh, bcs, A, rhs);
  
  if (rank == 0) { toc ("*** Boundary conditions ***"); }
  
  // Solve problem.
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  
  std::cout << "Solving linear system." << std::endl;
  
  mumps mumps_solver;
  
  if (rank == 0) { toc ("*** Solver - Initialize ***"); }
  
  // Set lhs.
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  
  std::vector<double> vals;
  std::vector<int> irow, jcol;
  
  A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
  
  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);
  
  if (rank == 0) { toc ("*** Solver - Set lhs ***"); }
  
  // Reduce rhs (so that rank 0 has the actual rhs).
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  
  std::vector<double> global_rhs(tmsh.num_global_nodes(), 0);
  MPI_Allreduce(rhs.data(), global_rhs.data(), rhs.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (rank == 0)
    mumps_solver.set_rhs (global_rhs);
  
  if (rank == 0) { toc ("*** Solver - Set rhs ***"); }
  
  // Solve.
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  mumps_solver.id.icntl[27] = 2; // Parallel analysis.
  mumps_solver.analyze ();
  if (rank == 0) { toc ("*** Solver - Analyze ***"); }
  
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  mumps_solver.factorize ();
  if (rank == 0) { toc ("*** Solver - Factorize ***"); }
  
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  mumps_solver.solve ();
  if (rank == 0) { toc ("*** Solver - Solve ***"); }
  
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  mumps_solver.cleanup ();
  if (rank == 0) { toc ("*** Solver - Cleanup ***"); }
  
  // Export solution.
  MPI_Barrier (MPI_COMM_WORLD); if (rank == 0) { tic (); }
  
  tmsh.vtk_export ("p4est_operator_segment_test");
  
  MPI_Bcast(global_rhs.data(), global_rhs.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  tmsh.octbin_export ("p4est_operator_segment_test_output", global_rhs);
  
  if (rank == 0) { toc ("*** Export ***"); }
  
  if (rank == 0) { print_timing_report (); }
  
  MPI_Finalize ();
  return 0;
}
