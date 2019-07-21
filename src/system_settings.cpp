#include "system_settings.h"


void
system_settings::read_param(std::string const path, std::string const namevar, int const num_param,
                              std::map<std::pair<int, int>, double> *output, GetPot fileData)
{
  for(int i = 0; i<num_param; i++)
  {
    for(int j = 0 ; j<num_param; j++)
    {
      std::string param = path+namevar;
      param +=std::to_string(i) + std::to_string(j);
      double value = fileData(param.c_str(),0.0);
      if (value != 0)
      {
        (*output)[std::make_pair(i,j)] = value;
      }
    }
  }
}

void
system_settings::read_param(std::string const path, std::string const namevar, int const num_param,
                              std::map<int, double> *output, GetPot fileData)
{
  for(int i = 0; i<num_param; i++)
  {
      std::string param = path+namevar;
      param +=std::to_string(i);
      double value = fileData(param.c_str(),0.0);
      if (value != 0)
      {
        (*output)[i] = value;
      }
  }
}

void
system_settings::read_param(std::string const path, std::string const namevar, int const num_param,
                              std::map<std::pair<int, int>, muparser_fun> *output, GetPot fileData)
{
  for(int i = 0; i<num_param; i++)
  {
    for(int j = 0 ; j<num_param; j++)
    {
      std::string param = path+namevar;
      param +=std::to_string(i) + std::to_string(j);
      std::string func = fileData(param.c_str(),"");
      if (func != "")
      {
        muparser_fun mu_fun (func,num_param);
        (*output)[std::make_pair(i,j)]=mu_fun;
      }
    }
  }
}

void
system_settings::read_param(std::string const path,std::string const namevar,int const num_param,
                              std::map<int, muparser_fun> *output, GetPot fileData)
{
  for(int i = 0; i<num_param; i++)
  {
      std::string param = path+namevar;
      param +=std::to_string(i);
      std::string func = fileData(param.c_str(),"");
      if (func != "")
      {
        muparser_fun mu_fun (func,num_param);
        (*output)[i]=mu_fun;
      }
  }
}


void
system_settings::read_from_file_reduced(std::string const filename){

  GetPot fileData(filename.c_str());
  muparser_fun zeroN ("0",N);

  ////////////////////////////READING FROM FILE/////////////////////////////////

  // function_estimator
  std::string param = "system/functions/fun_estimator/F";
  std::string func = fileData(param.c_str(),"x0");
  muparser_fun mu_fun (func,N);
  function_estimator=mu_fun;

  //Initial conditions
  for(int i = 0; i<N; i++)
  {
      std::string param = "system/functions/initial_cond/U";
      param +=std::to_string(i);
      std::string func = fileData(param.c_str(),"0");
      muparser_fun mu_fun (func,2);
      lambdas_init[i]=mu_fun;
    }

  // coeff_lap
  read_param("system/param/lapcoeff/","D",N,&coeff_lap,fileData);

  // coeff_time
  std::map<int, double> coeff_time;
  read_param("system/param/timecoeff/","A",N,&coeff_time,fileData);

  // coeff_adv
  std::map<std::pair<int, int>, double> coeff_adv;
  read_param("system/param/advcoeff/","MU",N,&coeff_adv,fileData);

  // lambdas_lik
  std::map<std::pair<int, int>, muparser_fun> lambdas_lik;
  read_param("system/functions/lin_f_angular/","L",N,&lambdas_lik,fileData);

  // lambdas_pik
  std::map<std::pair<int, int>, muparser_fun> lambdas_pik;
  read_param("system/functions/advection/","P",N,&lambdas_pik,fileData);

  // lambdas_bi
  std::map<int, muparser_fun> lambdas_bi;
  read_param("system/functions/lin_f_intercept/","B",N,&lambdas_bi,fileData);

  // fun rhs
  read_param("system/functions/forcing/","S",N,&fun_rhs,fileData);

  ////////////ASSEMBLING FUNCTIONS//////////////////////////////////////////////

  // lambdas_vec_reaz

  for( auto iter=coeff_time.begin(); iter!=coeff_time.end(); iter++){
    auto it=lambdas_lik.find(std::make_pair(iter->first,iter->first));
    if(it!=lambdas_lik.end()){
      muparser_fun mu_fun;
      mu_fun = it->second+(iter->second)/DELTAT;
      lambdas_vec_reaz[it->first]=mu_fun;
      lambdas_lik.erase(it);
    }
    else{
      muparser_fun mu_fun;
      mu_fun = zeroN + (iter->second)/DELTAT;
      lambdas_vec_reaz[std::make_pair(iter->first,iter->first)]=mu_fun;
    }
  }

  for (auto iter=lambdas_lik.begin(); iter!=lambdas_lik.end(); iter++)
    lambdas_vec_reaz[iter->first]=iter->second;


  // lambdas_vec_adv

  for (auto iter=coeff_adv.begin(); iter!=coeff_adv.end(); iter++){
    auto it= lambdas_pik.find(iter->first);
    if(it!=lambdas_pik.end()){
      muparser_fun mu_fun;
      mu_fun=(it->second)*(iter->second);
      lambdas_vec_adv[iter->first]=mu_fun;
    }
  }

  // lambdas_vec_forc

  for( auto iter=coeff_time.begin(); iter!=coeff_time.end(); iter++){
    std::string param ="x";
    param+=std::to_string(iter->first);
    muparser_fun funi (param,N);
    auto it=lambdas_bi.find(iter->first);
    if(it!=lambdas_bi.end()){
      muparser_fun mu_fun;
      mu_fun = funi*((iter->second)/DELTAT)-it->second;
      lambdas_vec_forc[it->first]=mu_fun;
      lambdas_bi.erase(it);
      coeff_time.erase(iter);
    }
    else{
      muparser_fun mu_fun;
      mu_fun = funi*((iter->second)/DELTAT);
      lambdas_vec_forc[iter->first]=mu_fun;
    }
  }

  for (auto iter=lambdas_bi.begin(); iter!=lambdas_bi.end(); iter++)
    lambdas_vec_forc[iter->first]=iter->second*(-1);

};

void
system_settings::read_from_file_complete(std::string const filename){

  GetPot fileData(filename.c_str());

  // Reading solver parameters

  N = fileData("system/solver_settings/NUM_EQUATIONS",4);
  NUM_ADAPT = fileData("system/solver_settings/NUM_ADAPT",1);
  NUM_NON_ADAPT = fileData("system/solver_settings/NUM_NON_ADAPT",1);
  REFINE_ITER = fileData("system/solver_settings/REFINE_ITER",1);
  TOL_EST = fileData("system/solver_settings/TOL_EST",1.);
  TOL_METRICS = fileData("system/solver_settings/TOL_METRICS",1.);
  PRINT_PROG = fileData("system/solver_settings/PRINT_PROGRESSIVE",1);
  PRINT_SOL = fileData("system/solver_settings/PRINT_SOL",1);
  PRINT_GRAD = fileData("system/solver_settings/PRINT_GRAD",1);
  PRINT_EST = fileData("system/solver_settings/PRINT_ESTIMATOR",1);
  DELTAT = fileData("system/solver_settings/DELTAT",1.);

  // Reading equations paramenters and functions
  read_from_file_reduced(filename);

};


void
system_settings::set_solver_parameter(int const N_,int const NUM_ADAPT_,int const NUM_NON_ADAPT_,int const REFINE_ITER_,
      double const DELTAT_,double const T_,int const PRINT_PROG_,int const PRINT_SOL_,int const PRINT_GRAD_,int const PRINT_EST_)
{
  N = N_;
  NUM_ADAPT = NUM_ADAPT_;
  NUM_NON_ADAPT = NUM_NON_ADAPT_;
  REFINE_ITER = REFINE_ITER_;
  PRINT_PROG = PRINT_PROG_;
  PRINT_SOL = PRINT_SOL_;
  PRINT_GRAD = PRINT_GRAD_;
  PRINT_EST = PRINT_EST_;
  DELTAT = DELTAT_;
};
