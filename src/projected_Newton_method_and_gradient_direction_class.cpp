
/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file projected_Newton_method_and_gradient_direction_class.cpp
  \brief interface for a nonlinear solver based on projected Newton method combined with projected gradient direction.
*/

#include <abstract_nonlinear_problem.h>
#include <projected_Newton_method_and_gradient_direction_class.h>
#include <bim_sparse.h>
#include <operators.h>
#include <nonlinear_solver.h>
#include <linear_solver.h>


void
projected_Newton_method_and_gradient_direction::set_problem
(abstract_nonlinear_problem *problem_)
{ problem = problem_; }

void
projected_Newton_method_and_gradient_direction::set_forcing_term
(abstract_forcing_term *forcing_)
{ forcing = forcing_; }

void
projected_Newton_method_and_gradient_direction::set_initial_guess
(std::vector<double> &initial_guess_)
{ initial_guess = &initial_guess_; }

void
projected_Newton_method_and_gradient_direction::set_bounds
(std::vector<double> &b1_ , std::vector<double> &b2_ )
{ b1 = b1_;
  b2 = b2_; }


void
projected_Newton_method_and_gradient_direction::projection
(std::vector<double> &x )
{
  for(int i = 0 ; i<x.size(); ++i)
    x[i] = std::max (std::min (x [i], b2 [i]), b1 [i]);

}


int
projected_Newton_method_and_gradient_direction::solve ()
{
  std::vector<int> ir, jc;
  std::vector<double> xa;
  std::vector<double> lin_initial_guess;
  sparse_matrix lhsT ;
  sparse_matrix mass_matrix;
  double theta_k = 1;
  std::vector<double> f_old, f_new, df_gap;

  #ifdef VERIFY_CONVERGENCE
    std::vector<double> linear_res, nonlinear_res, temp_res;
    std::vector<int> nonlinear_iter;
  #endif

  unsigned int n = b1.size(); 

  if (rank == 0)
    {
      problem->operator () (lhs, rhs, (*initial_guess));
      for (unsigned int i = 0; i < rhs.size (); ++i)
        residual_norm += rhs[i] * rhs[i];
      
      residual_norm = sqrt (residual_norm);
      std::cout<<"Initial residual "<< residual_norm<<std::endl;
      #ifdef VERIFY_CONVERGENCE
        nonlinear_res.push_back (residual_norm);
        linear_res.push_back (residual_norm);
        nonlinear_iter.push_back (0);
      #endif
       
    }
  
  MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (residual_norm < min_residual)
    {
      std::cout << "The Initial guess satisfies minimum residual."
                << std::endl;
      return 1;
    }
  
  if (rank == 0 && verbose >= 1)
    {
      // fout.open (filename.c_str ());
      
      std::cout << "Result of Non Linear Test "
                << "\nwill be written in "
                << filename
                << std::endl << std::endl;
   
      lin_initial_guess.assign (n, 0.0);
      
      if (lin_solver->solver_type () == "iterative")
	{
	  std::string type_of_iterative_method;
	  lin_solver->get_iterative_method(type_of_iterative_method);
	  std::cout <<"Linear solver used : ";
	  std::cout <<type_of_iterative_method << std::endl; 
	  std::cout <<std::endl;
	}
  
      printf ("%10.10s | \t%10.10s | \t%10.10s | \t%s | \t%10.10s\n",
	      "  Iterates", "    lambda",
	      "     ||F||", "       eta",
	      "      flag");

      printf ("__________________________________________________________________________\n");
    }
  
  iteration = 0; // DA AGGIUSTARE!!!
  bool FLAG = 0;
  do
    {
      ++iteration;
      if (rank ==0)
	{
             if (iteration > 1)
	       {
		 df_gap = lhs * rhs;
		 forcing_value = (*forcing)
		   (f_old, f_new, df_gap, forcing_value);
	       }
	     
             f_old.assign (rhs.size (), 0.0);
             for (unsigned int i = 0; i < rhs.size (); ++i)
	       f_old[i] = - rhs[i];

             lhs.aij (xa, ir, jc, lin_solver->get_index_base ());
             lin_solver->set_lhs_structure (lhs.rows (), ir, jc);
	}

      lin_solver->analyze ();
      
      if (rank == 0)
	lin_solver->set_lhs_data (xa);
      
      lin_solver->factorize ();
	  
      if (rank == 0)
	lin_solver->set_rhs (rhs);

      if (lin_solver->solver_type () == "iterative")
	{
	  lin_solver->set_tolerance (forcing_value);
	  lin_solver->set_initial_guess (lin_initial_guess);
	}

      lin_solver->solve ();   ///*///*///*///*///*///*///*///*///*///*
   
      if (rank == 0)
	{
	  std::vector<double> unew (initial_guess->size ());
	  for (unsigned int i = 0; i < rhs.size (); ++i)
	    unew[i] = rhs[i] + (*initial_guess)[i];
	  projection(unew);
	  double f_old_norm = residual_norm;
	  double f_new_norm = 0.0;
	  
	  (*problem) (f_new, unew);
	  
	  for (unsigned int i = 0; i < f_new.size (); ++i)
	    f_new_norm += f_new[i] * f_new[i];

	  f_new_norm = sqrt (f_new_norm);
	  
          #ifdef VERIFY_CONVERGENCE
	  temp_res.clear ();
	  temp_res = lhs * rhs;
	  for (unsigned int i = 0; i < temp_res.size (); ++i)
	    temp_res[i] += f_old[i];
	  
	  double temp_norm = 0.0;
	  for (unsigned int i = 0; i < temp_res.size (); ++i)
	    temp_norm += temp_res[i] * temp_res[i];

	  temp_norm = sqrt (temp_norm);

	  linear_res.push_back (temp_norm);
	  nonlinear_res.push_back (f_new_norm);
	  nonlinear_iter.push_back (iteration);
          #endif

	  unsigned int m = 0;
	  theta_k=1;  
	  while (f_new_norm >
		 (1 - t * theta_k*(1 - forcing_value)) * f_old_norm && m < max_back_it)
	    {
	      ///////// Theta Choice //////////////
	      std::vector<double> temp;
              double temp_norm = 0.0;

              temp = lhs * rhs;
              for (unsigned int i = 0; i < rhs.size (); ++i)
                temp_norm += f_old[i] * temp[i];

              double a = f_new_norm * f_new_norm -
                f_old_norm * f_old_norm - 2 * temp_norm;

              double b = 2* temp_norm;
              double c = f_old_norm * f_old_norm;

              theta_choice (a, b, c);
	     theta_k = theta ; 
              /////////

	     // theta_k=theta_k*theta;
         

                  for (unsigned int i = 0; i < rhs.size (); ++i)
		    {
		      rhs[i] *= theta;
		      unew[i] = rhs[i] + (*initial_guess)[i];
		    }
                  projection(unew);

                  (*problem) (f_new, unew);
		  f_new_norm = 0;  // l'ho aggiuto io !
                  for (unsigned int i = 0; i < f_new.size (); ++i)
		    f_new_norm += f_new[i] * f_new[i];

                  f_new_norm = sqrt (f_new_norm);
                  ++m;
	    }
            
         #ifdef VERIFY_CONVERGENCE
	  temp_res.clear ();
	  temp_res = lhs * rhs;
	  for (unsigned int i = 0; i < temp_res.size (); ++i)
	    temp_res[i] += f_old[i];
	  
	  temp_norm = 0.0;
	  for (unsigned int i = 0; i < temp_res.size (); ++i)
	    temp_norm += temp_res[i] * temp_res[i];

	  temp_norm = sqrt (temp_norm);
	  
	  linear_res.push_back (temp_norm);
	  nonlinear_res.push_back (f_new_norm);
	  nonlinear_iter.push_back (iteration);
         #endif
              std::cout<< "m " << m <<std::endl;
	  if (m > (max_back_it-1)) {FLAG = 1;}
	  else
	    {
	      for (unsigned int i = 0; i < rhs.size (); ++i)
		(*initial_guess)[i] = rhs[i] + (*initial_guess)[i];
	      projection(*initial_guess);
                 
	      (*problem) ( lhs,rhs, (*initial_guess));

	      residual_norm = 0.0;

	      for (unsigned int i = 0; i < rhs.size (); ++i)
		residual_norm += rhs[i] * rhs[i];
	      residual_norm = sqrt (residual_norm);
                  
	      if (verbose >=1)
		{
		  printf ("%10.5d |\t%10.5g |\t%10.5g |\t%10.5g |\t",
			  iteration, theta_k, residual_norm, forcing_value);
		  printf ("%10.10s\n", "PN");
		}
	    }
            
	} // rank ==0

      MPI_Bcast (&step_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&forcing_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
        
      if(FLAG == 1)
	{
          if (rank == 0 )
            { 
              // dato che il pattern di J^T non cambia, è meglio andare a modificare i valori,
	      // che crearla di nuovo ogni volta, quindi questo è da aggiustare 
              lhsT.reset ();
              lhsT.resize (n);

              if (lin_solver->solver_type () == "iterative")
		for (int i = 0 ; i <xa.size(); ++i)
		  lhsT[jc[i]][ir[i]] = xa[i] ;
	        
              else
		for (int i = 0 ; i <xa.size(); ++i)
		  lhsT[jc[i]-1][ir[i]-1] = xa[i] ;
		    
		
              rhs = lhsT*f_old;
	      for (int i = 0; i < rhs.size(); ++i)
		rhs[i] = - rhs[i];
	      std::vector<double> rhs_0 = rhs; 
              std::vector<double> unew (initial_guess->size ());
              for (unsigned int i = 0; i < rhs.size (); ++i)
                unew[i] = rhs[i] + (*initial_guess)[i];
	      projection(unew);
              double f_old_norm = residual_norm;
              double f_new_norm = 0.0;
	      
              (*problem) (f_new, unew);
	      
              for (unsigned int i = 0; i < f_new.size (); ++i)
		f_new_norm += f_new[i] * f_new[i];
	      
              f_new_norm = sqrt (f_new_norm);
              unsigned int m = 0;
              double val=0;
              std::vector<double> diff;
              diff.resize(n);

              for (int i=0 ; i <  n ; ++i)
		diff[i] = unew[i] - (*initial_guess)[i];
		
	      
              for (int i = 0 ; i< unew.size(); ++i)
		val += - sigma*rhs[i]*diff[i];

	      theta_k=1;
              while (0.5*f_new_norm*f_new_norm >
		     0.5*f_old_norm*f_old_norm + val && m < max_back_it )
		{
		  theta_k = theta_k * thetaPG;
		  
		  for (unsigned int i = 0; i < rhs.size (); ++i)
		    {
		      rhs[i] *= thetaPG;
		      unew[i] = rhs[i] + (*initial_guess)[i];
		    }
		  projection(unew);
		  
		  (*problem) (f_new, unew);
		  f_new_norm = 0;  // l'ho aggiuto io !
		  for (unsigned int i = 0; i < f_new.size (); ++i)
                    f_new_norm += f_new[i] * f_new[i];
		  
		  f_new_norm = sqrt (f_new_norm);
		  val=0;
		  
		  for (int i=0 ; i <  n ; ++i)
		    diff[i] = unew[i] - (*initial_guess)[i];
		    
		  for (int i = 0 ; i< unew.size(); ++i)
		    val += - sigma*rhs_0[i]*diff[i];
		   
		  ++m;
		}
	      /*
	      std :: cout << "new_norm = " <<0.5*f_new_norm*f_new_norm<<std::endl;
	      std :: cout << "old_norm = " <<0.5*f_old_norm*f_old_norm<<std::endl;
	      double min_unew = 300000; 
	      for (int i = 0 ; i < unew.size(); ++i)
		min_unew = std::min (min_unew, unew[i]);
	      std::cout << "min_new = " <<min_unew <<std::endl;
	      */
              FLAG=0;
              df_gap = lhs * rhs;
	      
	      for (unsigned int i = 0; i < rhs.size (); ++i)
		(*initial_guess)[i] = rhs[i] + (*initial_guess)[i];
	      projection(*initial_guess);
	      (*problem) ( lhs,rhs, (*initial_guess));

	      residual_norm = 0.0;

	      for (unsigned int i = 0; i < rhs.size (); ++i)
		residual_norm += rhs[i] * rhs[i];


	      residual_norm = sqrt (residual_norm);
              if (rank==0 && verbose >= 1)
                {
	           printf ("%10.5d |\t%10.5g |\t%10.5g |\t%10.5g |\t",
		      iteration, theta_k, residual_norm, forcing_value);
	           printf ("%10.10s\n", "PG"); 
	        }

            }// rank == 0
	  MPI_Bcast (&step_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Bcast (&residual_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Bcast (&forcing_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } // FLAG == 1
    }
  while (iteration <= max_iter
         && residual_norm > min_residual);

  if (rank == 0 && verbose >=1)
    // fout.close ();
  
  if (iteration > max_iter)
    {
      if (rank == 0)
        std::cout << "Maximum number of iterations reached."
                  << std::endl
                  << "Nonlinear solver not converged."
                  << std::endl;
     // return 0;
    }
  else
    {
      if (rank == 0)
        {
          std::cout << "\nNonlinear solver converged."
                    << std::endl;

          #ifdef VERIFY_CONVERGENCE
            std::cout << "Convergence of linear residual:" << std::endl;
            for (int i = 0; i < linear_res.size (); ++i)
              std::cout << nonlinear_iter[i] << ": "
                        << linear_res[i] << std::endl;

            std::cout << "Convergence of nonlinear residual:" << std::endl;
            for (int i = 0; i < nonlinear_res.size (); ++i)
              std::cout << nonlinear_iter[i] << ": "
                        << nonlinear_res[i] << std::endl;
          #endif
        }
      return 1;
    }
}

void
projected_Newton_method_and_gradient_direction::theta_choice
(double a, double b, double c)
{
  theta = - b / (2 * a);
  if (theta > theta_max || theta < theta_min || a < 0)
    {
      double f_min, f_max;
      f_min = (a * theta_min + b) * theta_min + c;
      f_max = (a * theta_max + b) * theta_max + c;
      if (f_min < f_max)
        theta = theta_min;
      else
        theta = theta_max;
    }
}

void
projected_Newton_method_and_gradient_direction::set_max_iterations_of_linear_solver
(int max_iteration)
{ lin_solver->set_max_iterations (max_iteration); }

void
projected_Newton_method_and_gradient_direction::set_initial_guess_of_linear_solver
(std::vector<double> &initial_guess)
{ lin_solver->set_initial_guess (initial_guess); }

void
projected_Newton_method_and_gradient_direction::set_initial_tolerance_of_linear_solver
(double initial_tolerance)
{
  lin_solver->set_tolerance (initial_tolerance);
  forcing_value = initial_tolerance;
}

void
projected_Newton_method_and_gradient_direction::set_iterative_method_of_linear_solver
(const std::string &iterative_method)
{
  lin_solver->set_iterative_method (iterative_method);
}

void
projected_Newton_method_and_gradient_direction::set_restart_iterations_of_linear_solver
(int restart_iterations)
{
  lin_solver->set_restart_iterations (restart_iterations);
}

void
projected_Newton_method_and_gradient_direction::set_options_iterative_method_of_linear_solver
(const std::string &options_iterative_method)
{
  lin_solver->set_options_iterative_method(options_iterative_method);
}

void
projected_Newton_method_and_gradient_direction::set_preconditioner_of_linear_solver
(const std::string &preconditioner)
{
  lin_solver->set_preconditioner (preconditioner);
}

void
projected_Newton_method_and_gradient_direction::set_convergence_condition_of_linear_solver
(const std::string &convergence_condition)
{
  lin_solver->set_convergence_condition (convergence_condition);
}

void
projected_Newton_method_and_gradient_direction::get_result_residual_norm
(double &residual_norm_)
{ residual_norm_ =  residual_norm; }

void
projected_Newton_method_and_gradient_direction::get_result_solution
(std::vector<double> &solution)
{ solution = *initial_guess; }

void
projected_Newton_method_and_gradient_direction::get_result_iterations
(int &iterations_)
{ iterations_ = iteration; }

void
projected_Newton_method_and_gradient_direction::cleanup ()
{ lin_solver->cleanup (); }

void
projected_Newton_method_and_gradient_direction::set_output_filename
(const std::string &filename_)
{ filename = filename_; }
