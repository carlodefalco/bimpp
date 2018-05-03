/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file lis_distributed_class.h
  \brief Interface for linear solver built with lis library.
*/

#ifndef HAVE_LIS_DISTRIBUTED_CLASS
#define HAVE_LIS_DISTRIBUTED_CLASS 1

#include <lis.h>
#include <string>
#include <bim_sparse.h>
#include <linear_solver.h>


/// Interface for linear solver built with lis library.
class lis_distributed: public linear_solver
{
private :

  LIS_SOLVER solver;
  LIS_MATRIX A;
  LIS_VECTOR b, x;
  LIS_INT iter;
  double time;
  bool initialized;
  
  LIS_INT *row,  *col;
  LIS_SCALAR *value;

  /// Stores rows of matrix in AIJ format.
  std::vector<int> irow;

  /// Stores columns of matrix in AIJ format
  std::vector<int> jcol;

  /// Pointer to data of matrix.
  double *data;

  /// Pointer to values of rhs.
  double *rhs;

  /// Stores values of ordering map.
  std::vector<int> ordering_map;

  /// Pointer to values of initial guess.
  double *initial_guess;

  /// It's true if solver have initial guess passed by user.
  bool have_initial_guess;

  /// \brief Maximum number of iterations.
  /// \details [default = 1000].
  int max_iter;

  /// \brief Number of iterations after which the iterative method is restarted.
  /// \details [default = 40].
  int restart_iterations;

  /// \brief Tolerance of linear solver.
  /// \details [default = 1e-12].
  double tolerance;

  /// \brief Name of preconditioner.
  /// \details [default = none].
  std::string preconditioner;
  
  /// \brief Type of convergence condition.
  /// \details [default = norm 2 of residual]
  std::string convergence_condition;

  std::vector<int> map_i_s, map_row_s;
  std::vector<int> map_n, map_nnz;

  int i_s, row_s;
  int n, nnz, n_row;

  /// Index base used by specific linear solver.
  static const int index_base = 0;

  int rank, size;

  /// Private helper functions
  
  void
  set_lhs_structure_aij
  (int n, std::vector<int> &ir, std::vector<int> &jc);

  int
  init_lis_objects ();

  int
  assemble_lis_matrix ();

  int
  invoke_lis_solver ();

  void
  cleanup_master ();

  void
  cleanup_slaves ();

  void
  destroy_lis_objects ();
  
public :

  /// \brief Name of iterative method.
  /// \details [default = bicg].
  std::string iterative_method;

  /// \brief Option string.
  std::string option_string;
  std::string options_iterative_method;
  bool option_string_set;

  bool verbose;
  
  /// Default costructor.
  lis_distributed (LIS_INT argc = 0, char * argv[] = NULL) :
    linear_solver ("LIS", "iterative"),
    initialized (false),
    row (0),
    col (0),
    value (0),
    data (0),
    rhs (0),
    initial_guess (0),
    have_initial_guess (false),
    max_iter (1000),
    restart_iterations(40),
    tolerance (1.0e-12),
    iterative_method ("bicg"),
    preconditioner ("none"),
    convergence_condition ("nrm2_r"),
    options_iterative_method (" "),
    option_string (""),
    option_string_set (false),
    verbose (true)
  {
    lis_initialize (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
  };

  /// Set-up the matrix structure.
  void
  set_lhs_structure
  (int n,
   std::vector<int> &ir,
   std::vector<int> &jc,
   matrix_format_t f)
  {
    if (f == aij)
      set_lhs_structure_aij (n, ir, jc);
    else
       std::cout << "Error: set_lhs_structure works only if the imput matrix has aij structure"
                <<std::endl;
     
  };

  
  /// Set matrix entries.
  void
  set_lhs_data (std::vector<double> &xa);

  /// Set the rhs.
  void
  set_rhs (std::vector<double> &rhs_)
  { rhs = &*rhs_.begin (); }

  /// Set the initial guess.
  void
  set_initial_guess (std::vector<double> &initial_guess_)
  {
    have_initial_guess = true;
    initial_guess = &*initial_guess_.begin ();
  }

  /// Prepare the solver.
  int
  factorize ()
  {
    init_lis_objects ();
    assert (assemble_lis_matrix () == 1);
    
    return 1;
  }

  /// Solve the system.
  int
  solve ()
  {

  MPI_Bcast (&have_initial_guess, 1, MPI_INT, 0, MPI_COMM_WORLD);
/*
  {
    MPI_Scatterv (&rhs[0], &map_n[0], &map_row_s[0], MPI_DOUBLE,
                  MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (have_initial_guess)  // da aggiustare!!!
      MPI_Scatterv (&initial_guess[0], &map_n[0], &map_row_s[0],
                    MPI_DOUBLE, MPI_IN_PLACE, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
*/
   int retval = invoke_lis_solver ();
 
   if (verbose && rank == 0)
     std::cout << std::endl
              << "Number of iterations = " << iter
              << std::endl
              << "Elapsed time = " << time << std::endl;
  
   MPI_Reduce (&rhs[0], MPI_IN_PLACE, n_row, MPI_DOUBLE, MPI_SUM, 0,
			      MPI_COMM_WORLD);

   return retval;


  }

  /// Cleanup memory.
  void
  cleanup ()
  {
    if (rank == 0)
      cleanup_master ();
    else
      cleanup_slaves ();

    destroy_lis_objects ();
    //lis_finalize ();
  }

  /// Set maximum number of iterations (default = 1000).
  void
  set_max_iterations (int max_iter_)
  {
    max_iter = max_iter_;
    option_string_set = false;
  }
  
  /// Set number of restart iterations (default = 40).
  void
  set_restart_iterations (int restart_iterations_)
  {
    restart_iterations = restart_iterations_;
    option_string_set = false;
  }

  /// Get maximum number of iterations.
  void
  get_max_iterations (int &max_iter_)
  {
    max_iter_ = max_iter;
    option_string_set = false;
  }

  /// Set tolerance of iterative method (default = 1.0e-12).
  void
  set_tolerance (double tol)
  {
    tolerance = tol;
    option_string_set = false;
  }

  /// Get tolerance of linear solver.
  void
  get_tolerance (double &tol)
  { tol = tolerance; }

  /// Get type of iterative method.
  void
  get_iterative_method (std::string &s)
  { s = iterative_method; }

  /// Get type of preconditioner.
  void
  get_preconditioner (std::string &s)
  { s = preconditioner; }

  /// get convergence condition of iterative method
  void
  get_convergence_condition (std::string &s)
  { s = convergence_condition; }

  /// Set type of iterative method (default = cg).
  void
  set_iterative_method (const std::string &s);

  /// Set options of iterative method .
  void
  set_options_iterative_method
  (const std::string &options_iterative_method_)
  {
    options_iterative_method = options_iterative_method_;
    option_string_set = false;
  }

  /// Set type of preconditioner (default = none).
  void
  set_preconditioner (const std::string &s);

  /// Set convergence condition of iterative method
  /// (default = nrm2_r)
  void
  set_convergence_condition (const std::string &s);

  /// LIS uses 0-based indexing
  inline int
  get_index_base ()
  { return index_base; }

};
#endif
