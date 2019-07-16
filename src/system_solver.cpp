#include "system_solver.h"

void
system_solver::set_init_cond( distributed_vector* init)
{
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
       quadrant != tmsh.end_quadrant_sweep ();
       ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){
            double xx=quadrant->p(0,ii);
            double yy=quadrant->p(1,ii);
            for(int i=0; i<N; i++){
              (*init)[ord[i](quadrant->gt (ii))] = lambdas_init[i](xx,yy);
            }
          }
          else
            {
              for(int i=0; i<N; i++){
              (*init)[ord[i](quadrant->gparent(0,ii))] +=0.;
              (*init)[ord[i](quadrant->gparent(1,ii))] +=0.;
            }

            }
        }
    }
    for(int i=0; i<N-1; i++)
      bim2a_solution_with_ghosts (tmsh, *init, replace_op, ord[i], false);
    bim2a_solution_with_ghosts (tmsh, *init, replace_op, ord[N-1]);
    TOC ("Set initial condition");
};

void
system_solver::print_sol_progressive( distributed_vector* to_print)
{
  char filename[255]="";
  for(int i=0; i<N; i++){
    sprintf(filename, "cahn_hilliard_u_%2.2d_%4.4d",i, count_progressive);
    tmsh.octbin_export (filename, *to_print, ord[i]);
  }
  count_progressive++;
};

void
system_solver::print_sol_analysis( distributed_vector* to_print, int m)
  {
    char filename[255]="";
    for(int i=0; i<N; i++){
      sprintf(filename, "cahn_hilliard_u_%2.2d_%4.4d_%2.2d",i, time_count,m);
      tmsh.octbin_export (filename, *to_print, ord[i]);
    }
  };

void
system_solver::print_grad_analysis( std::vector<gradient<distributed_vector>>* gradients, int m)
{
  char filename[255]="";
  for(int i=0; i<N; i++){
    sprintf(filename, "du_%2.2d_x_%4.4d_%2.2d",i,time_count,m);
    tmsh.octbin_export (filename, (*gradients)[i].first);
    sprintf(filename, "du_%2.2d_y_%4.4d_%2.2d",i,time_count,m);
    tmsh.octbin_export (filename, (*gradients)[i].second);
  }

};

void
system_solver::print_estim_analysis( std::vector<double>* estim_vec, int m)
{
  char filename[255]="";
  sprintf(filename, "estimator_%4.4d_%2.2d",time_count,m);
  tmsh.octbin_export_quadrant (filename, *estim_vec);
};

void
system_solver::obtain_global( distributed_vector* local, std::vector<distributed_vector>* global)
{
  for(int i=0; i<N; i++){
    std::vector<double> vec(gn_nodes);
    distributed_vector temp(gn_nodes);
    temp.get_owned_data ().assign (temp.get_owned_data ().size (), 0.0);
    for(int idx=0; idx<gn_nodes; ++idx)
      if(ord[i](idx)>=(*local).get_range_start () && ord[i](idx)<(*local).get_range_end ())
        vec[idx]=(*local)(ord[i](idx));
    MPI_Allreduce(MPI_IN_PLACE, vec.data(), gn_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int idx=0; idx<gn_nodes; ++idx)
      temp(idx)=vec[idx];
    (*global).push_back(temp);
}
};

void
system_solver::estimator_solution(std::vector<distributed_vector>* global, int m)
{
  std::vector<gradient<distributed_vector>> gradients;
  std::vector<q2_vec> u_stars;

  for(int i=0; i <N; i++){
    gradients.push_back(bim2c_quadtree_pde_recovered_gradient(tmsh, (*global)[i]));
    u_stars.push_back(bim2c_quadtree_pde_recovered_solution(tmsh, (*global)[i], gradients[i]));
  }
  if(PRINT_GRAD)
    print_grad_analysis(&gradients, m);
  std::vector<distributed_vector> temp=*global;

  auto estimator = [& u_stars, &temp] (tmesh::quadrant_iterator q)
     { return estimator_sol (q, u_stars[0], temp[0]); };
  tmsh.set_metrics_marker (estimator, TOL_EST, 4, 2, 2);

  // Compute estimator and h
  std::vector<double> estim_vec(ln_elements);

   double hx = 0, hy = 0,
   h = std::numeric_limits<double>::max (),
   global_h = 0;
   double est = 0, global_est = 0;

  for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep ();
        ++quadrant)
      {
        estim_vec[quadrant->get_forest_quad_idx ()] =
      estimator(quadrant)* std::sqrt (tmsh.num_global_quadrants ())
        / TOL_EST;

        hx = quadrant->p(0, 1) - quadrant->p(0, 0);
        hy = quadrant->p(1, 2) - quadrant->p(1, 0);

        h = std::min(h, std::sqrt(hx*hx + hy*hy));

        est += std::pow(estimator(quadrant), 2);
      }
    if(PRINT_EST)
      print_estim_analysis(&estim_vec,m);

  MPI_Reduce(&h, &global_h, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&est, &global_est, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  global_est = std::sqrt(global_est);

  nnodes.push_back (tmsh.num_global_nodes ());
  h_step.push_back (global_h);
  estim.push_back (global_est);
};

void
system_solver::interpolation(distributed_vector* old_vec, distributed_vector* new_vec)
{
  for (int i=0; i<N; i++){
    interpolate_vector (tmsh, *old_vec, *new_vec, ord[i]);
  }
  for(int i=0; i<N-1; i++)
    bim2a_solution_with_ghosts (tmsh, *new_vec, replace_op, ord[i], false);
  bim2a_solution_with_ghosts (tmsh, *new_vec, replace_op, ord[N-1]);
};

void
system_solver::linear_solution(distributed_sparse_matrix* A, distributed_vector* rhs,
                                mumps* lin_solver, std::vector<double>* xa, std::vector<int>* ir, std::vector<int>* jc)
{




  // Matrix update
  TIC ();
  (*A).aij_update (*xa, *ir, *jc, lin_solver->get_index_base ());
  lin_solver->set_distributed_lhs_data (*xa);
  TOC ("set LHS data");


  // Factorization
  TIC ();
  std::cout << "lin_solver->factorize () = " << lin_solver->factorize () << std::endl;
  TOC ("solver factorize");


  // Set RHS data
  TIC ();
  lin_solver->set_rhs_distributed (*rhs);
  TOC ("set RHS data");


  // Solution
  TIC ();
  std::cout << "lin_solver->solve () = " << lin_solver->solve () << std::endl;
  TOC ("solver solve");

};

void
system_solver::assemble_lin_sys (distributed_sparse_matrix* A, distributed_vector* rhs)
{
  TIC();
  (*A).assemble();
  (*rhs).assemble();
  TOC("Communicate A and b");
};

void
system_solver::mesh_value_update()
{
  gn_nodes = tmsh.num_global_nodes ();
  ln_nodes = tmsh.num_owned_nodes ();
  ln_elements = tmsh.num_local_quadrants ();
};

void
system_solver::mesh_init()
{
  // Generation of initial mesh
  constexpr p4est_topidx_t simple_conn_num_vertices = 4;
  constexpr p4est_topidx_t simple_conn_num_trees = 1;
  const double simple_conn_p[simple_conn_num_vertices*2] =
    {0., 0., 1., 0.,  1., 1., 0., 1.};
  const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] =
    {1, 2, 3, 4, 1};

  tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                            simple_conn_t, simple_conn_num_trees);
  int recursive = 1;
  tmsh.set_refine_marker (uniform_refinement);
  tmsh.refine (recursive);
};

void
system_solver::element_evaluate(std::vector<std::vector<double>>* cont_vec, std::vector<double>* coeff)
{
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep ();
      ++quadrant)
    {
      for(int i=0; i<N*N; i++)
      {
        (*cont_vec)[i][quadrant->get_forest_quad_idx ()]=(*coeff)[i];
      }
    }
  TOC ("compute constant coefficients and initial conditions");

};

void
system_solver::nodes_evaluate(std::vector<distributed_vector>* cont_vec,
                              std::vector<std::function<double(std::vector<double>)>>* lambdas_vec,
                              distributed_vector* uold)
{
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep ();
      ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){
            std::vector<double> u_nodes;
            for(int j=0; j<N; j++)
              u_nodes.push_back((*uold)[ord[j](quadrant->gt (ii))]);
            for(int i =0; i<cont_vec->size(); i++)
              (*cont_vec)[i][quadrant->gt (ii)] = (*lambdas_vec)[i](u_nodes);
          }
          else
          {
            for(int i=0; i<cont_vec->size(); i++){
              (*cont_vec)[i][quadrant->gparent(0,ii)] +=0;
              (*cont_vec)[i][quadrant->gparent(1,ii)] +=0;
            }
          }

       }
     }
  for(int i=0; i<cont_vec->size(); i++)
    (*cont_vec)[i].assemble(replace_op);
  TOC ("compute constant coefficients and initial conditions");
};

void
system_solver::assemble_matrix(distributed_sparse_matrix* A, distributed_vector* rhs, distributed_vector* uold)
{
  // IN qualche modo otteniamo coeff_lap,  lambdas_vec_adv, lambda_vec_reaz, lambda_vec_rhs, lambdas_vec_forc, fun_rhs


  std::vector<double> ones (ln_elements,1);


//////////////////////////////////BEGIN TEST///////////////////////////////////
  std::function<double(std::vector<double>)> zero;
  zero = [&] (std::vector<double> x) {return 0.0;};

  std::vector<double> coeff_lap(N*N,0.0);
  coeff_lap[1]=1;
  coeff_lap[4]=EPSU*EPSU;
  coeff_lap[11]=1;
  coeff_lap[15]=EPSV*EPSV;


  std::vector<std::function<double(std::vector<double>)>> lambdas_vec_adv;
  for(int i=0; i<N*N; i++){
    lambdas_vec_adv.push_back(zero);
  }

  std::vector<std::function<double(std::vector<double>)>> lambdas_vec_reaz;
  for(int i=0; i<N*N; i++){
    lambdas_vec_reaz.push_back(zero);
  }

  lambdas_vec_reaz[0]= [&] (std::vector<double> x) {return -TAUU/DELTAT;};
  lambdas_vec_reaz[4]= [&] (std::vector<double> x) {return 1.5*x[0]*x[0]-0.5;};
  lambdas_vec_reaz[6]= [&] (std::vector<double> x) {return ALPHA+BETA*x[2];};
  lambdas_vec_reaz[5]= [&] (std::vector<double> x) {return 1;};
  lambdas_vec_reaz[10]= [&] (std::vector<double> x) {return -SIGMA-TAUV/DELTAT;};
  lambdas_vec_reaz[14]= [&] (std::vector<double> x) {return 1.5*x[2]*x[2]+2*BETA*x[0]-0.5;};
  lambdas_vec_reaz[12]= [&] (std::vector<double> x) {return ALPHA;};
  lambdas_vec_reaz[15]= [&] (std::vector<double> x) {return 1;};

  std::vector<std::function<double(std::vector<double>)>> lambdas_vec_forc;
  for(int i=0; i<N; i++){
    lambdas_vec_forc.push_back(zero);
  }

  lambdas_vec_forc[0]= [&] (std::vector<double> x) {return -x[0]*TAUU/DELTAT;};
  lambdas_vec_forc[1]= [&] (std::vector<double> x) {return 0.5*x[0]*x[0]*x[0]+0.5*x[0];};
  lambdas_vec_forc[2]= [&] (std::vector<double> x) {return -x[2]*TAUV/DELTAT;};
  lambdas_vec_forc[3]= [&] (std::vector<double> x) {return 0.5*x[2]*x[2]*x[2]+0.5*x[2];};

  std::vector<std::function<double(double t, double x, double y)>> fun_rhs(N);
  fun_rhs[2] = [&] (double t, double x, double y) {return -VBAR*SIGMA;};
  fun_rhs[0] = [&] (double t, double x, double y) {return 0.0;};
  fun_rhs[1] = [&] (double t, double x, double y) {return 0.0;};
  fun_rhs[3] = [&] (double t, double x, double y) {return 0.0;};



////////////////////////////////////END TEST//////////////////////////////////

  std::vector<std::vector<double>> lap(N*N, std::vector<double> (ln_elements));
  element_evaluate(&lap, &coeff_lap);
  for( int i=0; i<N; i++)
    for( int k=0; k<N; k++)
      bim2a_laplacian(tmsh, lap[i*N+k], *A, ord[i], ord[k]);


  std::vector<distributed_vector> upwind(N*N, distributed_vector(ln_nodes));
  nodes_evaluate(&upwind, &lambdas_vec_adv, uold);
  for( int i=0; i<N; i++)
    for( int k=0; k<N; k++)
      bim2a_advection_upwind(tmsh, upwind[i*N+k], *A, ord[i], ord[k]);

  std::vector<distributed_vector> reaz(N*N, distributed_vector(ln_nodes));
  nodes_evaluate(&reaz, &lambdas_vec_reaz, uold);
  for( int i=0; i<N; i++)
    for( int k=0; k<N; k++)
       bim2a_reaction(tmsh, ones, reaz[i*N+k], *A, ord[i], ord[k]);


  std::vector<distributed_vector> forc(N, distributed_vector(ln_nodes));
  nodes_evaluate(&forc, &lambdas_vec_forc, uold);
    for( int i=0; i<N; i++)
          bim2a_rhs (tmsh, ones, forc[i], *rhs, ord[i]);

  std::vector<distributed_vector> forc_add(N, distributed_vector(ln_nodes));
  TIC ();
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
        quadrant != tmsh.end_quadrant_sweep ();
          ++quadrant)
        {
        for (int ii = 0; ii < 4; ++ii)
          {
            if (! quadrant->is_hanging (ii)){
                double xx= quadrant->p(0,ii);
                double yy= quadrant->p(1,ii);
                for(int i =0; i<N; i++)
                    forc_add[i][quadrant->gt (ii)] = fun_rhs[i](time_count*DELTAT, xx, yy);
                }
                else
                {
                for(int i=0; i<N; i++){
                    forc_add[i][quadrant->gparent(0,ii)] +=0;
                    forc_add[i][quadrant->gparent(1,ii)] +=0;
                  }
                }

              }
            }
      for(int i=0; i<N; i++)
        forc_add[i].assemble(replace_op);
    TOC ("Add forcing term to rhs");
    for( int i=0; i<N; i++)
          bim2a_rhs (tmsh, ones, forc_add[i], *rhs, ord[i]);
};

void
system_solver::solver_analysis(distributed_sparse_matrix* A, mumps * lin_solver, std::vector<double> *xa, std::vector<int> *ir, std::vector<int> *jc )
{
  lin_solver->set_lhs_distributed ();
  A->aij (*xa, *ir, *jc, lin_solver->get_index_base ());
  lin_solver->set_distributed_lhs_structure (A->rows (), *ir, *jc);
  std::cout << "lin_solver->analyze () = "<< lin_solver->analyze () << std::endl;
};

void
system_solver::obtain_solution (mumps* lin_solver, distributed_vector* new_sol)
{
  distributed_vector result = lin_solver->get_distributed_solution ();
  for (int idx = new_sol->get_range_start (); idx < new_sol->get_range_end (); ++idx)
    (*new_sol) (idx) = result (idx);
  for(int i=0; i<N-1; i++)
    bim2a_solution_with_ghosts (tmsh, *new_sol, replace_op, ord[i], false);
   bim2a_solution_with_ghosts (tmsh, *new_sol, replace_op, ord[N-1]);
};

void
system_solver::solve ()
{

  mesh_init();
  mesh_value_update();

  // Typedef
  using q1_vec  = q1_vec<distributed_vector>;

  // Time
  if(rank==0)
    std::cout<<"TIME= "<<0<<std::endl;

  // Initial data containers
  q1_vec init (ln_nodes * N);
  init.get_owned_data ().assign (init.get_owned_data ().size (), 0.0);

  // Set initial conditions
  set_init_cond(&init);

  // Save initial conditions on old mesh
  if(PRINT_PROG)
    print_sol_progressive(&init);
  if(PRINT_SOL)
    print_sol_analysis(&init,0);


  // Copy local u(0) in u(0)_g
  TIC();
  std::vector<q1_vec> global;
  obtain_global(&init, &global);
  TOC("Obtaining global solution");

  // Now I have u(0)_g on M(0)

  // Estimator for u(0)
  TIC();
  estimator_solution(&global,0);
  TOC("Compute estim and h")

  // Refine and obtain new parameters
  TIC();
  tmsh.metrics_refine (TOL_METRICS);
  mesh_value_update();
  TOC("refine");

  // Interpolate u(0) on M(1)
  TIC();
  q1_vec sold (ln_nodes * N);
  sold.get_owned_data ().assign (sold.get_owned_data ().size (), 0.0);
  interpolation(&init, &sold);
  TOC("Interpolation");

  // Now I have u(0) on M(1)

  // Print u(0) on M(1)
  if(PRINT_PROG)
    print_sol_progressive(&sold);
  if(PRINT_SOL)
    print_sol_analysis(&sold,1);


  // Copy local u(0) in u(0)_g
  TIC();
  global.clear();
  obtain_global(&sold, &global);
  TOC("Obtaining global solution");


  // Now I have u(0)_g on M(1)
  TIC();
  estimator_solution(&global,1);
  TOC("Compute estim and h")

   ///////////////////////////////TIME CYCLE///////////////////////////////////


   // Vectors to use linear solver
  std::vector<double> xa;
  std::vector<int> ir, jc;


  for (int adapt =0; adapt<NUM_ADAPT; adapt++){

    // Declare Matrix and solver for first NUM_NON_ADAPT steps
    TIC();
    mumps *lin_solver = new mumps ();

    distributed_sparse_matrix A;
    A.set_ranges (ln_nodes * N);

    q1_vec soldd (ln_nodes * N);
    soldd.get_owned_data ().assign (soldd.get_owned_data ().size (), 0.0);



    q1_vec rhs (ln_nodes * N);
    rhs.get_owned_data ().assign (rhs.get_owned_data ().size (), 0.0);

    xa.clear();
    ir.clear();
    jc.clear();


    assemble_matrix(&A, &rhs, &sold);
    assemble_lin_sys(&A, &rhs);

    // Solver analysis
    TIC ();
    solver_analysis(&A, lin_solver, &xa, &ir, &jc);
    TOC ("solver analysis");

    // NUM_NON_ADAPT time steps
    for (int j =0 ; j< NUM_NON_ADAPT; j++){
         time_count++;

        // Print curent time
        if(rank==0)
          std::cout<<"TIME= "<<DELTAT*time_count<<std::endl;


          // Reset containers
        TIC();
        A.reset ();
        rhs.get_owned_data ().assign (rhs.get_owned_data ().size (), 0.0);
        rhs.assemble (replace_op);
        TOC("Resetting")

        assemble_matrix(&A,&rhs,&sold);
        assemble_lin_sys(&A, &rhs);
        linear_solution(&A, &rhs, lin_solver, &xa, &ir, &jc);

        // Copy solution
        TIC();
        if(j==NUM_NON_ADAPT-1)
            soldd=sold;
        obtain_solution (lin_solver, &sold);
        TOC("Obtaining solution");

        // Print solution
        TIC();
        if(PRINT_PROG)
          print_sol_progressive(&sold);
        if(PRINT_SOL)
          print_sol_analysis(&sold,0);
        TOC("Exporting solution");
    }

    lin_solver->cleanup();

    // After some time steps on fixed mesh we exit from the cycle with
    // u(s) on M(r) and u(s-1) on M(r)

    // REFINE_ITER steps of refinement

    for(int m=0; m <REFINE_ITER; m++){

      // Obtaining global solution
      TIC();
      global.clear();
      obtain_global(&sold, &global);
      TOC("Obtaining global solution");

      TIC();
      estimator_solution(&global,m);
      TOC("Compute estim and h")
      TIC();

      // Refine
      TIC();
      tmsh.metrics_refine (TOL_METRICS);
      mesh_value_update();
      TOC("refine");

      // Interpolate u(s-1) on M(r+1)
      TIC();
      q1_vec soldd_interp (ln_nodes * N);
      soldd_interp.get_owned_data ().assign (soldd_interp.get_owned_data ().size (), 0.0);
      interpolation(&soldd, &soldd_interp);
      soldd=soldd_interp; // saving u(s-1) on M(r+1) in soldd for next refinement step
      TOC("Interpolation");

       // Now we have u(s-1) on M(r+1)

       // ONE TIME STEP to obtain u(s) on M(r+1), computed through u(s-1) on M(r+1)

       // Prepare containers to solve
       TIC();
       q1_vec rhs (ln_nodes * N);
       rhs.get_owned_data ().assign (rhs.get_owned_data ().size (), 0.0);


       mumps *lin_solver = new mumps ();

       xa.clear();
       ir.clear();
       jc.clear();

       distributed_sparse_matrix A;
       A.set_ranges (ln_nodes * N);
       assemble_matrix(&A, &rhs, &soldd_interp);
       assemble_lin_sys(&A,&rhs);

       TOC("Containers construction");

       // Solver analysis
       TIC ();
       solver_analysis(&A, lin_solver, &xa, &ir, &jc);
       TOC ("solver analysis");

       linear_solution(&A, &rhs, lin_solver, &xa, &ir, &jc);


       // Copy solution
       TIC();
       obtain_solution(lin_solver, &soldd_interp);
       TOC("Obtaining solution");


       // Save solution
       TIC();
       if(PRINT_PROG)
        print_sol_progressive(&soldd_interp);
       if(PRINT_SOL)
        print_sol_analysis(&soldd_interp, m+1);
       sold=soldd_interp; // Saving u(s) on M(r+1) in sold for next refinement step
       lin_solver->cleanup ();
       TOC("Exporting solution");

     }
  }


   // Print reports
   if (rank == 0)
     for (unsigned step = 0; step < nnodes.size(); ++step)
        std::cout << "Step " << step << ", #nodes: "
                  << nnodes[step] << ", h: "
                  << h_step[step] << ", estimator: "
                  << estim[step] << std::endl;

  MPI_Barrier (MPI_COMM_WORLD);
  if (rank == 0) { print_timing_report (); }

};
