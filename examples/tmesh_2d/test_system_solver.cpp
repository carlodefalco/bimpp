#include <system_solver.h>


int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  std::vector<std::function<double(double x, double y)>> initcond;
  initcond.push_back([&] (double x, double y) { return std::sin(10*x*y); });
  initcond.push_back([&] (double x, double y) { return 0.0; });
  initcond.push_back([&] (double x, double y) { return std::cos(10*(x-y))*x*y; });
  initcond.push_back([&] (double x, double y) { return 0.0; });

  std::vector<std::function<double(double x)>> lineariz;
  lineariz.push_back([&] (double x) { return -1.5*x*x+0.5; });
  lineariz.push_back([&] (double x) { return 0.5*x*(x*x+1);});



  system_solver sys(rank, initcond, lineariz, 4, 50, 4, 2, 1e-4, 1e5, 0.005, 5, 0.05, 0.05, 1, 100, 0.04, -0.9, 0.0114559, 100,1,0,0,0);
  sys.solve();

  MPI_Finalize ();

  return 0;
}
