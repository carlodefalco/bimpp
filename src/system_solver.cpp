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
            std::vector<double> temp = {xx,yy};
            for(auto iter=settings.lambdas_init.begin(); iter!=settings.lambdas_init.end(); iter++)
              (*init)[ord[iter->first](quadrant->gt (ii))] = (iter->second)(temp);
          }
          else
            {
              for(int i=0; i<settings.N; i++){
              (*init)[ord[i](quadrant->gparent(0,ii))] +=0.;
              (*init)[ord[i](quadrant->gparent(1,ii))] +=0.;
            }

            }
        }
    }
    for(int i=0; i<settings.N-1; i++)
      bim2a_solution_with_ghosts (tmsh, *init, replace_op, ord[i], false);
    bim2a_solution_with_ghosts (tmsh, *init, replace_op, ord[settings.N-1]);
    TOC ("Set initial condition");
};

void
system_solver::print_sol_progressive( distributed_vector* to_print)
{
  char filename[255]="";
  for(int i=0; i<settings.N; i++){
    sprintf(filename, "cahn_hilliard_u_%2.2d_%4.4d",i, count_progressive);
    tmsh.octbin_export (filename, *to_print, ord[i]);
  }
  count_progressive++;
};

void
system_solver::print_sol_analysis( distributed_vector* to_print, int m)
{
    char filename[255]="";
    for(int i=0; i<settings.N; i++){
      sprintf(filename, "cahn_hilliard_u_%2.2d_%4.4d_%2.2d",i, time_count,m);
      tmsh.octbin_export (filename, *to_print, ord[i]);
    }
  };

void
system_solver::print_grad_analysis( std::vector<gradient<distributed_vector>>* gradients, int m)
{
  char filename[255]="";
  for(int i=0; i<settings.N; i++){
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
  for(int i=0; i<settings.N; i++)
  {
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
system_solver::estimator_solution(std::vector<distributed_vector>* global,int m)
{
  std::vector<gradient<distributed_vector>> gradients;
  std::vector<q2_vec> u_stars;


  for(int i=0; i <settings.N; i++){
    gradients.push_back(bim2c_quadtree_pde_recovered_gradient(tmsh, (*global)[i]));
    u_stars.push_back(bim2c_quadtree_pde_recovered_solution(tmsh, (*global)[i], gradients[i]));
  }

  if(settings.PRINT_GRAD)
    print_grad_analysis(&gradients, m);

  std::vector<distributed_vector> temp=*global;
  muparser_fun mu_fun = settings.function_estimator;
  std::vector<int> index_sol_estim=mu_fun.get_var_index(0,settings.N);
  std::vector<int> index_grad_estim=mu_fun.get_var_index(settings.N,2*settings.N);

  auto estimator = [& u_stars, &temp, &mu_fun, &index_sol_estim, &index_grad_estim, &gradients] (tmesh::quadrant_iterator q)
  {
    std::vector<double> sol_estimators(temp.size(),0.0);
    for(int i=0; i<index_sol_estim.size(); i++)
      sol_estimators[index_sol_estim[i]]=(estimator_sol (q, u_stars[index_sol_estim[i]], temp[index_sol_estim[i]]));

    std::vector<double> grad_estimators(temp.size(),0.0);
    for(int i=0; i<index_grad_estim.size(); i++)
     grad_estimators[index_grad_estim[i]]=(estimator_grad(q, gradients[index_grad_estim[i]], temp[index_grad_estim[i]]));

    sol_estimators.insert(sol_estimators.end(), grad_estimators.begin(), grad_estimators.end());

    return mu_fun(sol_estimators);
  };

  tmsh.set_metrics_marker (estimator, settings.TOL_EST, 4, 2, 2);

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
        estim_vec[quadrant->get_forest_quad_idx ()] = estimator(quadrant)*
        std::sqrt (tmsh.num_global_quadrants ())
        / settings.TOL_EST;

        hx = quadrant->p(0, 1) - quadrant->p(0, 0);
        hy = quadrant->p(1, 2) - quadrant->p(1, 0);

        h = std::min(h, std::sqrt(hx*hx + hy*hy));

        est += std::pow(estimator(quadrant), 2);
  }

  if(settings.PRINT_EST)
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
  for (int i=0; i<settings.N; i++){
    interpolate_vector (tmsh, *old_vec, *new_vec, ord[i]);
  }
  for(int i=0; i<settings.N-1; i++)
    bim2a_solution_with_ghosts (tmsh, *new_vec, replace_op, ord[i], false);
  bim2a_solution_with_ghosts (tmsh, *new_vec, replace_op, ord[settings.N-1]);
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
  (*A).assemble();
  (*rhs).assemble();
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
system_solver::element_evaluate(std::vector<double>* cont_vec, double value)
{
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep ();
      ++quadrant)
    {
      (*cont_vec)[quadrant->get_forest_quad_idx ()]=value;
      }

};

void
system_solver::nodes_evaluate(distributed_vector* cont_vec,
                              muparser_fun* lambdas_vec,
                              distributed_vector* uold)
{
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
      quadrant != tmsh.end_quadrant_sweep ();
      ++quadrant)
    {
      for (int ii = 0; ii < 4; ++ii)
        {
          if (! quadrant->is_hanging (ii)){
            std::vector<double> u_nodes;
            for(int j=0; j<settings.N; j++)
              u_nodes.push_back((*uold)[ord[j](quadrant->gt (ii))]);
            (*cont_vec)[quadrant->gt (ii)] = (*lambdas_vec)(u_nodes);
          }
          else
          {

              (*cont_vec)[quadrant->gparent(0,ii)] +=0;
              (*cont_vec)[quadrant->gparent(1,ii)] +=0;
            }
          }

       }

  (*cont_vec).assemble(replace_op);
};

void
system_solver::assemble_matrix(distributed_sparse_matrix* A, distributed_vector* rhs, distributed_vector* uold)
{
  std::vector<double> ones (ln_elements,1);
  TIC();
  for(auto iter=settings.coeff_lap.begin(); iter!=settings.coeff_lap.end(); iter++){
    std::vector<double> lap(ln_elements);
    element_evaluate(&lap, iter->second);
    bim2a_laplacian(tmsh, lap, *A, ord[(iter->first).first], ord[(iter->first).second]);
  }

  for(auto iter=settings.lambdas_vec_reaz.begin(); iter!=settings.lambdas_vec_reaz.end(); iter++){
    distributed_vector reaz(ln_nodes);
    nodes_evaluate(&reaz, &(iter->second), uold);
    bim2a_reaction(tmsh, ones, reaz, *A, ord[(iter->first).first], ord[(iter->first).second]);
  }

  for(auto iter=settings.lambdas_vec_adv.begin(); iter!=settings.lambdas_vec_adv.end(); iter++){
    distributed_vector adv(ln_nodes);
    nodes_evaluate(&adv, &(iter->second), uold);
    bim2a_advection_upwind(tmsh, adv, *A, ord[(iter->first).first], ord[(iter->first).second]);
  }

  for(auto iter=settings.lambdas_vec_forc.begin(); iter!=settings.lambdas_vec_forc.end(); iter++){
    distributed_vector forc(ln_nodes);
    nodes_evaluate(&forc, &(iter->second), uold);
    bim2a_rhs(tmsh, ones, forc, *rhs, ord[iter->first]);
  }
  TOC("computing coefficients")

  TIC ();
  for (auto iter=settings.fun_rhs.begin(); iter!=settings.fun_rhs.end(); iter++)
  {
    distributed_vector forc_add (ln_nodes);
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
          quadrant != tmsh.end_quadrant_sweep ();
            ++quadrant)
    {
          for (int ii = 0; ii < 4; ++ii)
          {
            if (! quadrant->is_hanging (ii)){
              double xx= quadrant->p(0,ii);
              double yy= quadrant->p(1,ii);
              std::vector<double> temp= {time_count*settings.DELTAT, xx, yy};
              forc_add[quadrant->gt (ii)] = (iter->second)(temp);
            }
            else
              forc_add[quadrant->gparent(0,ii)] +=0;
          }

      }
    forc_add.assemble(replace_op);
    bim2a_rhs (tmsh, ones, forc_add, *rhs, ord[iter->first]);
  }
  TOC ("add forcing term to rhs");
};

void
system_solver::fake_assemble_matrix(distributed_sparse_matrix* A, distributed_vector* rhs)
{
  TIC();
  std::vector<double> ones_elements (ln_elements,1.0);

  distributed_vector ones_nodes(ln_nodes);
  ones_nodes.get_owned_data().assign(ones_nodes.get_owned_data().size(),0.0);
  for (auto quadrant = tmsh.begin_quadrant_sweep ();
        quadrant != tmsh.end_quadrant_sweep ();
          ++quadrant)
        {
        for (int ii = 0; ii < 4; ++ii)
          {
            if (! quadrant->is_hanging (ii)){
                    ones_nodes[quadrant->gt (ii)] = 1.0;
                }
            else{
                ones_nodes[quadrant->gparent(0,ii)] +=0;
                ones_nodes[quadrant->gparent(1,ii)] +=0;
              }

          }

        }
  ones_nodes.assemble(replace_op);
  TOC("assemble ones");

  TIC();
  for(auto iter=settings.coeff_lap.begin(); iter!=settings.coeff_lap.end(); iter++)
    bim2a_laplacian(tmsh, ones_elements, *A, ord[(iter->first).first], ord[(iter->first).second]);


  for(auto iter=settings.lambdas_vec_reaz.begin(); iter!=settings.lambdas_vec_reaz.end(); iter++)
    bim2a_reaction(tmsh, ones_elements, ones_nodes , *A, ord[(iter->first).first], ord[(iter->first).second]);


  for(auto iter=settings.lambdas_vec_adv.begin(); iter!=settings.lambdas_vec_adv.end(); iter++)
    bim2a_advection_upwind(tmsh, ones_nodes, *A, ord[(iter->first).first], ord[(iter->first).second]);

  for(auto iter=settings.lambdas_vec_forc.begin(); iter!=settings.lambdas_vec_forc.end(); iter++)
    bim2a_rhs(tmsh, ones_elements, ones_nodes, *rhs, ord[iter->first]);

  for (auto iter=settings.fun_rhs.begin(); iter!=settings.fun_rhs.end(); iter++)
    bim2a_rhs(tmsh, ones_elements, ones_nodes, *rhs, ord[iter->first]);

  TOC("fake assmeble");
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
  for(int i=0; i<settings.N-1; i++)
    bim2a_solution_with_ghosts (tmsh, *new_sol, replace_op, ord[i], false);
   bim2a_solution_with_ghosts (tmsh, *new_sol, replace_op, ord[settings.N-1]);
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
  q1_vec init (ln_nodes * settings.N);
  init.get_owned_data ().assign (init.get_owned_data ().size (), 0.0);

  // Set initial conditions
  set_init_cond(&init);

  // Save initial conditions on old mesh
  if(settings.PRINT_PROG)
    print_sol_progressive(&init);
  if(settings.PRINT_SOL)
    print_sol_analysis(&init,0);


  // Copy local u(0) in u(0)_g
  TIC();
  std::vector<q1_vec> global;
  obtain_global(&init, &global);
  TOC("obtaining global solution");

  // Now I have u(0)_g on M(0)

  // Estimator for u(0)
  TIC();
  estimator_solution(&global,0);
  TOC("compute estim and h")

  // Refine and obtain new parameters
  TIC();
  tmsh.metrics_refine (settings.TOL_METRICS);
  mesh_value_update();
  TOC("refine");

  // Interpolate u(0) on M(1)
  TIC();
  q1_vec sold (ln_nodes * settings.N);
  sold.get_owned_data ().assign (sold.get_owned_data ().size (), 0.0);
  interpolation(&init, &sold);
  TOC("interpolation");

  // Now I have u(0) on M(1)

  // Print u(0) on M(1)
  if(settings.PRINT_PROG)
    print_sol_progressive(&sold);
  if(settings.PRINT_SOL)
    print_sol_analysis(&sold,1);


  // Copy local u(0) in u(0)_g
  TIC();
  global.clear();
  obtain_global(&sold, &global);
  TOC("obtaining global solution");


  // Now I have u(0)_g on M(1)
  TIC();
  estimator_solution(&global,1);
  TOC("compute estim and h")

   ///////////////////////////////TIME CYCLE///////////////////////////////////


   // Vectors to use linear solver
  std::vector<double> xa;
  std::vector<int> ir, jc;


  for (int adapt =0; adapt<settings.NUM_ADAPT; adapt++){

    // Declare Matrix and solver for first NUM_NON_ADAPT steps
    TIC();
    mumps *lin_solver = new mumps ();

    distributed_sparse_matrix A;
    A.set_ranges (ln_nodes * settings.N);

    q1_vec soldd (ln_nodes * settings.N);
    soldd.get_owned_data ().assign (soldd.get_owned_data ().size (), 0.0);

    q1_vec rhs (ln_nodes * settings.N);
    rhs.get_owned_data ().assign (rhs.get_owned_data ().size (), 0.0);

    xa.clear();
    ir.clear();
    jc.clear();

    fake_assemble_matrix(&A, &rhs);

    TIC();
    assemble_lin_sys(&A, &rhs);
    TOC("communicate A and b");
    // Solver analysis
    TIC ();
    solver_analysis(&A, lin_solver, &xa, &ir, &jc);
    TOC ("solver analysis");

    // NUM_NON_ADAPT time steps
    for (int j =0 ; j< settings.NUM_NON_ADAPT; j++)
    {
         time_count++;

        // Print curent time
        if(rank==0)
          std::cout<<"TIME= "<<settings.DELTAT*time_count<<std::endl;


          // Reset containers
        TIC();
        A.reset ();
        rhs.get_owned_data ().assign (rhs.get_owned_data ().size (), 0.0);
        rhs.assemble (replace_op);
        TOC("Resetting")

        assemble_matrix(&A,&rhs,&sold);

        TIC();
        assemble_lin_sys(&A, &rhs);
        TOC("communicate A and b");

        linear_solution(&A, &rhs, lin_solver, &xa, &ir, &jc);

        // Copy solution
        TIC();
        if(j==settings.NUM_NON_ADAPT-1)
            soldd=sold;
        obtain_solution (lin_solver, &sold);
        TOC("Obtaining solution");

        // Print solution
        TIC();
        if(settings.PRINT_PROG)
          print_sol_progressive(&sold);
        if(settings.PRINT_SOL)
          print_sol_analysis(&sold,0);
        TOC("Exporting solution");
    }

    lin_solver->cleanup();

    // After some time steps on fixed mesh we exit from the cycle with
    // u(s) on M(r) and u(s-1) on M(r)

    // REFINE_ITER steps of refinement

    for(int m=0; m <settings.REFINE_ITER; m++){

      // Obtaining global solution
      TIC();
      global.clear();
      obtain_global(&sold, &global);
      TOC("obtaining global solution");

      TIC();
      estimator_solution(&global,m);
      TOC("compute estim and h")

      // Refine
      TIC();
      tmsh.metrics_refine (settings.TOL_METRICS);
      mesh_value_update();
      TOC("refine");

      // Interpolate u(s-1) on M(r+1)
      TIC();
      q1_vec soldd_interp (ln_nodes * settings.N);
      soldd_interp.get_owned_data ().assign (soldd_interp.get_owned_data ().size (), 0.0);
      interpolation(&soldd, &soldd_interp);
      soldd=soldd_interp; // saving u(s-1) on M(r+1) in soldd for next refinement step
      TOC("interpolation");

       // Now we have u(s-1) on M(r+1)

       // ONE TIME STEP to obtain u(s) on M(r+1), computed through u(s-1) on M(r+1)

       // Prepare containers to solve
       TIC();
       q1_vec rhs (ln_nodes * settings.N);
       rhs.get_owned_data ().assign (rhs.get_owned_data ().size (), 0.0);


       mumps *lin_solver = new mumps ();

       xa.clear();
       ir.clear();
       jc.clear();

       distributed_sparse_matrix A;
       A.set_ranges (ln_nodes * settings.N);
       assemble_matrix(&A, &rhs, &soldd_interp);

       TIC();
       assemble_lin_sys(&A,&rhs);
       TOC("communicate A and b");

       TOC("containers construction");

       // Solver analysis
       TIC ();
       solver_analysis(&A, lin_solver, &xa, &ir, &jc);
       TOC ("solver analysis");

       linear_solution(&A, &rhs, lin_solver, &xa, &ir, &jc);


       // Copy solution
       TIC();
       obtain_solution(lin_solver, &soldd_interp);
       TOC("obtaining solution");


       // Save solution
       TIC();
       if(settings.PRINT_PROG)
        print_sol_progressive(&soldd_interp);
       if(settings.PRINT_SOL)
        print_sol_analysis(&soldd_interp, m+1);
       sold=soldd_interp; // Saving u(s) on M(r+1) in sold for next refinement step
       lin_solver->cleanup ();
       TOC("exporting solution");
     }

     TIC();
     global.clear();
     obtain_global(&sold, &global);
     TOC("obtaining global solution");

     TIC();
     estimator_solution(&global,settings.REFINE_ITER);
     TOC("compute estim and h")
     TIC();
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
