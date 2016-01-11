
#ifndef RRE_H
#define RRE_H

#include "bim_config.h"
#include <cassert>
#include <vector>
#include <iostream>

#define DGELSY_F77 F77_FUNC(dgelsy, DGELSY)
#define DCOPY_F77 F77_FUNC(dcopy, DCOPY)
#define DAXPY_F77 F77_FUNC(daxpy, DAXPY)
#define DGEMV_F77 F77_FUNC(dgemv, DGEMV)
#define DNRM2_F77 F77_FUNC(dnrm2, DNRM2)
extern "C"
{
  void
  DGELSY_F77(const int *M, const int *N, const int *NRHS, double *A,
             const int *LDA, double *B, int *LDB, int *JPVT,
             const double *RCOND, int *RANK, double *WORK, 
             const int *LWORK, int *INFO);

  void
  DCOPY_F77(const int* n,
            const double *dx,
            const int* incx,
            double *dy,
            const int* incy);

  void
  DAXPY_F77(const int* n,
            const double *da,
            const double *dx,
            const int* incx,
            double *dy,
            const int* incy);

  void
  DGEMV_F77 (const char* trans,
             const int* m,
             const int* n,
             const double* alpha,
             const double* a,
             const int* lda,
             const double* x,
             const int* incx,
             const double* beta,
             double *y,
             const int* incy);

  double
  DNRM2_F77 (const int* m,             
             const double* x,
             const int* incx);
  
}


class rre
{
private:

  int ninit;
  int nskip;
  int rank;
  int m;
  int nu;
  int nv;
  const int nrhs = 1;
  double rcond = 1.0e-17;
  int iter = 0;
  int lwork = -1;
  int info;
  int erank;
  
  double *X = 0;
  double *U = 0;
  double *V = 0;
  double *work = 0;
  int    *jpvt = 0;

  std::vector<double> nrm;
  
  int
  perform_extrapolation (double* x);
  
public:

  enum stage_t { init, skip, store } stage = init;

  rre (int size_, int ninit_ = 10,
       int nskip_ = 0, int rank_ = 5);

  ~rre ()
  {
    delete [] X;
    delete [] U;
    delete [] V;
    delete [] work;
  }
  
  int
  extrapolate (std::vector<double> &x);

  stage_t
  get_stage () { return stage; }

  int
  get_iter () { return iter; }

  void
  print_X ()
  {
    std::cout << "X = [" << std::endl;
    for (int jj = 0; jj < m; ++jj)
      {
        for (int ii = 0; ii < rank; ++ii)
          std::cout << X[(ii * m) + jj] << ", ";
        std::cout << std::endl;
      }
    std::cout << "];" << std::endl;
  }

  void
  print_U ()
  {
    std::cout << "U = [" << std::endl;
    for (int jj = 0; jj < m; ++jj)
      {
        for (int ii = 0; ii < rank-1; ++ii)
          std::cout << U[(ii * m) + jj] << ", ";
        std::cout << std::endl;
      }
    std::cout << "];" << std::endl;
  }

  void
  print_V ()
  {
    std::cout << "V = [" << std::endl;
    for (int jj = 0; jj < m; ++jj)
      {
        for (int ii = 0; ii < rank-2; ++ii)
          std::cout << V[(ii * m) + jj] << ", ";
        std::cout << std::endl;
      }
    std::cout << "];" << std::endl;
  }

  void
  print_nrm ()
  {
    std::cout << "nrm = [" << std::endl;
    for (unsigned int jj = 0; jj < nrm.size (); ++jj)
      std::cout << nrm[jj] << ";" << std::endl;
    std::cout << "];" << std::endl;
  }
  
};

#endif
