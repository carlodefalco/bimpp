#include "bim_config.h"
#include "rre_aux.h"
#include "rre.h"

const char trans = 'n';
const int step = 1;
const double mone = -1.0;
const double pone = +1.0;

rre::rre (int size_, int ninit_,
          int nskip_, int rank_)
  : ninit (ninit_),
    nskip (nskip_), rank (rank_),
    m (size_), nu (rank_ - 1),
    nv (rank_ - 2)
{
  assert (rank >= 3);
  double dtmp = 0.0;    
    
  X = new double [m * rank];
  U = new double [m * nu];
  V = new double [m * nv];
  jpvt = new int [nv] ();
      
  DGELSY_F77 (&m, &nv, &nrhs, V,   
              &m, U, &m, jpvt,
              &rcond, &erank, &dtmp,  
              &lwork, &info);

  /*
  std::cout << "lwork = " << lwork << std::endl;
  std::cout << "work[0] = " << dtmp << std::endl;
  std::cout << "info = " << info << std::endl;
  */

  lwork = static_cast<int> (dtmp);
  work = new double [lwork];
  
};

int
rre::extrapolate (std::vector<double> &x)
{
  
  assert (static_cast<unsigned int>(m) == x.size ());
  int info = 0;
  
  switch (stage)
    {
    case init :
      if (! (iter < ninit))
        {
          iter = -1;
          stage = store;
        }
      break;
      
    case store :

      DCOPY_F77 (&m, &x[0], &step,
                 &X[iter * m], &step);
      if (iter > 0)
        {              
          DCOPY_F77 (&m, &X[iter * m], &step,                         
                     &U[(iter-1) * m], &step);
          
          DAXPY_F77 (&m, &mone, &X[(iter-1) * m], &step,
                     &U[(iter-1) * m], &step);

          nrm.push_back (DNRM2_F77 (&m, &U[(iter-1) * m], &step));
            
          if (iter > 1)
            {
              DCOPY_F77 (&m, &U[(iter-1) * m], &step,
                         &V[(iter-2) * m], &step);
              
              DAXPY_F77 (&m, &mone, &U[(iter-2) * m], &step,
                         &V[(iter-2) * m], &step);
            }
        }
      if (iter == rank-1)
        {

          info = perform_extrapolation (&x[0]);
          //std::cout << "info = " << info << std::endl;
          
          iter = -1;
          if (nskip > 0)
            stage = skip; 
        }

      break;

    case skip :
      if (! (iter < nskip))
        {
          iter = -1;
          stage = store;
        }
      break;
    default :
      std::cerr << "impossible state reached!!!" << std::endl;
    }

  iter++;
  return info;
  
};


int
rre::perform_extrapolation (double *x)
{

  // U will be overwritten, but
  // we will need to use U(:, 1)
  // later. As U(:, nu) will not
  // be needed after extrapolation,
  // we store U(:, 1) into U(:, nu)
  // and multiply pinv(V) * U(:, nu)
  DCOPY_F77 (&m, U, &step,                         
             &U[m * (nu-1)], &step);
  
  DGELSY_F77 (&m, &nv, &nrhs, V,   
              &m, &U[m * (nu-1)], &m, jpvt,
              &rcond, &erank, work,  
              &lwork, &info);
  
  if (info < 0)
    return info;

  // here we need to use only rank-2
  // columns of U in the product.
  // The last column which would be unused,
  // contains the result of the least-squares
  // system above.
  DGEMV_F77 (&trans, &m, &nv, &mone,
             U, &m, &U[m * (nu-1)], &step, &pone,
             X, &step);
  
  DCOPY_F77 (&m, X, &step,                         
             x, &step);  
    
  return 0;
};
