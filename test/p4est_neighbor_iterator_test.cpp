#include <iostream>
#include <vector>

#include <bim_timing.h>
#include <tmesh.h>

#include <simple_connectivity_2d.h>


static int
uniform_refinement (tmesh::quadrant_iterator quadrant)
{ return 1; }

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

  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);

  recursive = 0; partforcoarsen = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive, partforcoarsen);
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive, partforcoarsen);
  
  for (int iproc = 0; iproc < size; ++iproc)
    {
      if ((tmsh.num_local_quadrants () > 0) && (iproc == rank))
        {
          for (auto quadrant = tmsh.begin_quadrant_sweep ();
               quadrant != tmsh.end_quadrant_sweep ();
               ++quadrant)
            {
              for (auto neighbor = quadrant->begin_neighbor_sweep ();
                   neighbor != quadrant->end_neighbor_sweep ();
                   ++neighbor)
                {
                  std::cout
                    << "Element "
                    << quadrant->get_global_quad_idx ()
                    << ", vertices: ";

                  for (int node = 0; node < 4; ++node)
                    {
                      if (!quadrant->is_hanging(node))
                        std::cout << quadrant->gt(node) << ", ";
                      else
                        std::cout
                          << "(" << quadrant->gparent (0, node)
                          << ", " << quadrant->gparent (1, node)
                          << "), ";
                    }          
                    
                  std::cout
                    << "neighbor "
                    << neighbor->get_global_quad_idx ()
                    << " (face " << neighbor.get_face_idx ()
                    << ")"
                    << ", vertices: ";

                  for (int node = 0; node < 4; ++node)
                    {
                      if (!neighbor->is_hanging(node))
                        std::cout << neighbor->gt(node) << ", ";
                      else
                        std::cout
                          << "(" << neighbor->gparent (0, node)
                          << ", " << neighbor->gparent (1, node)
                          << "), ";
                    }
                  std::cout << std::endl;
                }
              std::cout << std::endl;
            }
        }
      MPI_Barrier (mpicomm);
    }
  
  tmsh.vtk_export ("p4est_neighbor_iterator_test");

  MPI_Finalize ();
  return 0;

}
