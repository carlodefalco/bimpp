#include <system_solver.h>


int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);

  std::string file = "cahn_hilliard.txt";
  //std::string file = "cahn_hilliard_system.txt";

  // Reading from file example
  system_solver sys(file);
  sys.solve();

  // Reading a system_setting object example
  //system_settings problem(4,2,2,1,1e-4,1e5,0.005,5,1,0,0,0,file);
  //system_solver sys;
  //system_solver.set_problem(problem);


  // Reading partially from file
  //system_solver sys(4,6,2,1,1e-4,1e5,0.005,5,1,0,0,0,file);
  //sys.solve();

  MPI_Finalize ();

  return 0;
}
