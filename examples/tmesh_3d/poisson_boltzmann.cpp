#include <mpi.h>

#include "pb_class.h"
#include "pqr_parser.cpp"

static char filename[255];

int
main (int argc, char **argv)
{

  MPI_Init (&argc, &argv);
  
  int                   recursive, partforcoarsen, balance;
  MPI_Comm              mpicomm = MPI_COMM_WORLD;  
  int                   rank, size;
  tmesh_3d              tmsh;

  MPI_Comm_rank (mpicomm, &rank);
  MPI_Comm_size (mpicomm, &size);
 
  poisson_boltzmann pb;
  
  if (pb.parse_options (argc, argv))
    return 1;

  std::ifstream inputfile (pb.pqrfilename);
  read_atoms_from_pqr (inputfile, pb.atoms);
  inputfile.close ();
  
  if (rank == 0)
  {
     std::cout << "Atom : " << std::endl;
     write_atoms_to_pqr (std::cout, pb.atoms);
     pb.print_options ();
  }

  TIC ();
  if (pb.mesh_shape == 1)
     pb.create_mesh ();
  else if (pb.mesh_shape == 0)
     pb.create_cubic_mesh ();
  else 
  {
     std::cerr << "Invalid mesh shape selected" << std::endl;
     return 1;
  }
  TOC ("create_mesh");

  TIC ();
  pb.init_tmesh ();
  TOC ("init_tmesh");
  
  TIC ();
  pb.refine_surface ();
  TOC ("refine surface");

  TIC ();
  pb.create_markers ();
  TOC ("create element markers");
  
  TIC ();
  pb.export_ls_tmesh ();
  TOC ("export ls tmesh");

  TIC ();
  pb.export_marked_tmesh ();
  TOC ("export marked tmesh");
  
  TIC ();
  if (pb.linear_solver_name == "mumps")
     pb.mumps_compute_electric_potential ();
  else if (pb.linear_solver_name == "lis")
     pb.lis_compute_electric_potential ();
  else 
  {
     std::cerr << "Invalid linear solver selected" << std::endl;
     return 1;
  }
  TOC ("compute electric potential");
  
  TIC ();
  pb.export_p4est ();
  TOC ("export p4est");
    
  if (rank == 0) { print_timing_report(); }
  MPI_Barrier (mpicomm);
  
  MPI_Finalize ();
  return 0;
  
}

