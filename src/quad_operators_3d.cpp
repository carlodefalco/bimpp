#include "quad_operators_3d.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <set>
#include <limits>
#include <iomanip>

double
hm (const double & a, const double & b)
{
  return 2. / (1. / a + 1. / b);
}

// MPI_User_function.
static void replace(double *invec, double *inoutvec,
                    int *len, MPI_Datatype *dtype)
{
  for (int i = 0; i < *len; ++i)
    if (invec[i] != 0 && inoutvec[i] == 0)
      inoutvec[i] = invec[i];
}

// 4-points Gauss quadature nodes and weights (in [0, 1]).
static constexpr double gn[4] =
  {6.94318442029737e-02, 3.30009478207572e-01,
   6.69990521792428e-01, 9.30568155797026e-01};
static constexpr double gw[4] =
  {1.73927422568727e-01, 3.26072577431273e-01,
   3.26072577431273e-01, 1.73927422568727e-01};

// Transform nodes from [0, 1] to [x[0], x[1]].
static inline double
xformx (const double *x, const double X)
{ return (X * (x[1] - x[0]) + x[0]); }

// Transform weights from [0, 1] to [x[0], x[1]].
static inline double
xformw (const double *x, const double w)
{ return (w * (x[1] - x[0])); }

