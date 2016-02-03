#ifndef RRE_AUX_H
#define RRE_AUX_H

#include "bim_config.h"
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

#endif
