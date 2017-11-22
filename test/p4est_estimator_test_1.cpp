#include <bim_sparse.h>
#include <mumps_class.h>
#include <quad_operators.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>

#include <vector>
#include <cassert>

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
    write_example_connectivity ("p4est_estimator_test_1.octbin.gz");

  tmsh.read_connectivity ("p4est_estimator_test_1.octbin.gz");
  
  tmsh.set_refine_marker (uniform_refinement);
  recursive = 0; partforcoarsen = 1;
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.refine (recursive, partforcoarsen);
  
  tmsh.vtk_export ("p4est_estimator_test_1");
  
  // Assemble advection-diffusion matrix.
  sparse_matrix A;
  A.resize(tmsh.num_global_nodes());
  
  double lambda = 25;
  std::vector<double> alpha(tmsh.num_local_quadrants (), 1);
  std::vector<double> psi(tmsh.num_local_nodes (), 0);
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          psi[quadrant->t(ii)] = lambda * (quadrant->p(0, ii) + quadrant->p(1, ii));
        }
    }
  
  bim2a_advection_diffusion (tmsh, alpha, psi, A);
  
  // Assemble right-hand side.
  std::vector<double> rhs(tmsh.num_global_nodes (), 0);
  
  std::vector<double> f(tmsh.num_local_quadrants (), 0);
  std::vector<double> g(tmsh.num_local_nodes (), 0);
  
  bim2a_rhs (tmsh, f, g, rhs);
  
  // Set boundary conditions.
  func u_ex =
    [lambda] (double x, double y)
    { return (exp(lambda * x) - 1) / (exp(lambda) - 1) *
             (exp(lambda * y) - 1) / (exp(lambda) - 1); };
             
  dirichlet_bcs bcs;
  for (int i = 0; i < 4; ++i)
    bcs.push_back (std::make_tuple(0, i, u_ex));
  
  bim2a_dirichlet_bc (tmsh, bcs, A, rhs);
  
  // Solve problem.
  std::cout << "Solving linear system.";
  
  mumps mumps_solver;
  
  std::vector<double> vals;
  std::vector<int> irow, jcol;
  
  A.aij(vals, irow, jcol, mumps_solver.get_index_base ());
  
  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (A.rows (), irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);
  
  // Reduce rhs (so that rank 0 has the actual rhs).
  std::vector<double> global_rhs(tmsh.num_global_nodes(), 0);
  MPI_Allreduce(rhs.data(), global_rhs.data(), rhs.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (rank == 0)
    mumps_solver.set_rhs (global_rhs);
  
  // Solve.
  mumps_solver.analyze ();
  mumps_solver.factorize ();
  mumps_solver.solve ();
  mumps_solver.cleanup ();
  
  // Export solution.
  std::vector<double> u(global_rhs);
  MPI_Bcast(u.data(), u.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  tmsh.octbin_export ("p4est_estimator_test_1_output", u);
  
  std::cout << " Done." << std::endl;
  
  // Compute reconstructed gradient.
  std::cout << "Computing reconstructed gradient.";
  
  double hx = 0;
  double hy = 0;
  
  std::vector<double> du_x_star(tmsh.num_global_nodes(), 0);
  std::vector<double> du_y_star(tmsh.num_global_nodes(), 0);
  
  int node_n = 0, node_side = 0;
  std::vector<double> du_x, weights_x;
  std::vector<double> du_y, weights_y;
  
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      hx = quadrant->p(0, 1) - quadrant->p(0, 0);
      hy = quadrant->p(1, 2) - quadrant->p(1, 0);
      
      // Loop over vertices of current quadrant that
      // are non-hanging and have not been processed yet.
      for (int node = 0; node < 4; ++node)
        {
          if (quadrant->is_hanging(node)
              || du_x_star[quadrant->t(node)] != 0
              || du_y_star[quadrant->t(node)] != 0)
            break;
          
          du_x.clear(); weights_x.clear();
          du_y.clear(); weights_y.clear();
          
          switch (node)
            {
              case 0:
                du_x.push_back((u[quadrant->t(1)] - u[quadrant->t(0)]) / hx);
                du_y.push_back((u[quadrant->t(2)] - u[quadrant->t(0)]) / hy);
                break;
              case 1:
                du_x.push_back((u[quadrant->t(1)] - u[quadrant->t(0)]) / hx);
                du_y.push_back((u[quadrant->t(3)] - u[quadrant->t(1)]) / hy);
                break;
              case 2:
                du_x.push_back((u[quadrant->t(3)] - u[quadrant->t(2)]) / hx);
                du_y.push_back((u[quadrant->t(2)] - u[quadrant->t(0)]) / hy);
                break;
              case 3:
                du_x.push_back((u[quadrant->t(3)] - u[quadrant->t(2)]) / hx);
                du_y.push_back((u[quadrant->t(3)] - u[quadrant->t(1)]) / hy);
                break;
            }
          
          weights_x.push_back(1 / hx);
          weights_y.push_back(1 / hy);
           
          // Loop over face neighbors of current quadrant.
          for (auto neighbor = quadrant->begin_neighbor_sweep ();
               neighbor != quadrant->end_neighbor_sweep ();
               ++neighbor)
            {
              // Skip missing neighbors.
              if (neighbor->get_global_quad_idx() == quadrant->get_global_quad_idx())
                continue;
              
              // Check if neighbor contains current vertex ("node").
              for (node_n = 0; node_n < 4; ++node_n)
                if (neighbor->t(node_n) == quadrant->t(node))
                  break;
              
              // If not, switch to the next neighbor.
              if (node_n == 4)
                continue;
              
              hx = neighbor->p(0, 1) - neighbor->p(0, 0);
              hy = neighbor->p(1, 2) - neighbor->p(1, 0);
              
              switch (node_n)
                {
                  case 0:
                    if (node == 1)
                      du_x.push_back((u[neighbor->t(1)] - u[neighbor->t(0)]) / hx);
                    if (node == 2)
                      du_y.push_back((u[neighbor->t(2)] - u[neighbor->t(0)]) / hy);
                    break;
                  case 1:
                    if (node == 0)
                      du_x.push_back((u[neighbor->t(1)] - u[neighbor->t(0)]) / hx);
                    if (node == 3)
                      du_y.push_back((u[neighbor->t(3)] - u[neighbor->t(1)]) / hy);
                    break;
                  case 2:
                    if (node == 3)
                      du_x.push_back((u[neighbor->t(3)] - u[neighbor->t(2)]) / hx);
                    if (node == 0)
                      du_y.push_back((u[neighbor->t(2)] - u[neighbor->t(0)]) / hy);
                    break;
                  case 3:
                    if (node == 2)
                      du_x.push_back((u[neighbor->t(3)] - u[neighbor->t(2)]) / hx);
                    if (node == 1)
                      du_y.push_back((u[neighbor->t(3)] - u[neighbor->t(1)]) / hy);
                    break;
                }
              
              if (weights_x.size() < du_x.size())
                weights_x.push_back(1 / hx);
              
              if (weights_y.size() < du_y.size())
                weights_y.push_back(1 / hy);
            }
          
          // If on any vertical boundary/interface.
          if (du_x.size() < 2)
            {
              for (auto neighbor = quadrant->begin_neighbor_sweep ();
                   neighbor != quadrant->end_neighbor_sweep ();
                   ++neighbor)
                {
                  // Skip missing neighbors.
                  if (neighbor->get_global_quad_idx() == quadrant->get_global_quad_idx())
                    continue;
                  
                  // Check if neighbor contains the opposite vertex
                  // of the horizontal side containing "node".
                  switch (node)
                    {
                      case 0:
                        node_side = 1;
                        break;
                      case 1:
                        node_side = 0;
                        break;
                      case 2:
                        node_side = 3;
                        break;
                      case 3:
                        node_side = 2;
                        break;
                    }
                  
                  for (node_n = 0; node_n < 4; ++node_n)
                    if (neighbor->t(node_n) == quadrant->t(node_side))
                      break;
                  
                  // If not, switch to the next neighbor.
                  if (node_n == 4)
                    continue;
                  
                  hx = neighbor->p(0, 1) - neighbor->p(0, 0);
                  
                  switch (node_n)
                    {
                      case 0:
                        if (node_side == 1)
                          du_x.push_back((u[neighbor->t(1)] - u[neighbor->t(0)]) / hx);
                        break;
                      case 1:
                        if (node_side == 0)
                          du_x.push_back((u[neighbor->t(1)] - u[neighbor->t(0)]) / hx);
                        break;
                      case 2:
                        if (node_side == 3)
                          du_x.push_back((u[neighbor->t(3)] - u[neighbor->t(2)]) / hx);
                        break;
                      case 3:
                        if (node_side == 2)
                          du_x.push_back((u[neighbor->t(3)] - u[neighbor->t(2)]) / hx);
                        break;
                    }
                  
                  if (weights_x.size() < du_x.size())
                    {
                      weights_x[0] += 2 / hx;
                      weights_x.push_back(-1 / hx);
                    }
                }
            }
          
          // If on any horizontal boundary/interface.
          if (du_y.size() < 2)
            {
              for (auto neighbor = quadrant->begin_neighbor_sweep ();
                   neighbor != quadrant->end_neighbor_sweep ();
                   ++neighbor)
                {
                  // Skip missing neighbors.
                  if (neighbor->get_global_quad_idx() == quadrant->get_global_quad_idx())
                    continue;
                  
                  // Check if neighbor contains the opposite vertex
                  // of the vertical side containing "node".
                  switch (node)
                    {
                      case 0:
                        node_side = 2;
                        break;
                      case 1:
                        node_side = 3;
                        break;
                      case 2:
                        node_side = 0;
                        break;
                      case 3:
                        node_side = 1;
                        break;
                    }
                  
                  for (node_n = 0; node_n < 4; ++node_n)
                    if (neighbor->t(node_n) == quadrant->t(node_side))
                      break;
                  
                  // If not, switch to the next neighbor.
                  if (node_n == 4)
                    continue;
                  
                  hy = neighbor->p(1, 2) - neighbor->p(1, 0);
                  
                  switch (node_n)
                    {
                      case 0:
                        if (node_side == 2)
                          du_y.push_back((u[neighbor->t(2)] - u[neighbor->t(0)]) / hy);
                        break;
                      case 1:
                        if (node_side == 3)
                          du_y.push_back((u[neighbor->t(3)] - u[neighbor->t(1)]) / hy);
                        break;
                      case 2:
                        if (node_side == 0)
                          du_y.push_back((u[neighbor->t(2)] - u[neighbor->t(0)]) / hy);
                        break;
                      case 3:
                        if (node_side == 1)
                          du_y.push_back((u[neighbor->t(3)] - u[neighbor->t(1)]) / hy);
                        break;
                    }
                  
                  if (weights_y.size() < du_y.size())
                    {
                      weights_y[0] += 2 / hy;
                      weights_y.push_back(-1 / hy);
                    }
                }
            }
          
          for (unsigned int ix = 0; ix < du_x.size(); ++ix)
            {
              du_x_star[quadrant->t(node)] += du_x[ix] * weights_x[ix];
            }
          
          du_x_star[quadrant->t(node)] /= std::accumulate(weights_x.begin(), weights_x.end(), 0.0);
          
          for (unsigned int iy = 0; iy < du_y.size(); ++iy)
            {
              du_y_star[quadrant->t(node)] += du_y[iy] * weights_y[iy];
            }
          
          du_y_star[quadrant->t(node)] /= std::accumulate(weights_y.begin(), weights_y.end(), 0.0);
        }
    }
  
  tmsh.octbin_export ("p4est_estimator_test_1_output_du_x", du_x_star);
  tmsh.octbin_export ("p4est_estimator_test_1_output_du_y", du_y_star);
  
  std::cout << " Done." << std::endl;
  
  MPI_Finalize ();
  
  return 0;
}
