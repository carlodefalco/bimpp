#include <cmath>
#include <limits>
#include <mosfet_doping_2d.h>

static double
gaussian (const double x, const double c, const double s)
{ return (exp (- std::pow (x - c, 2) / (2 * std::pow (s, 2)))); }

static inline double
bool1 (const double x, const double L)
{ return (((x > L / 3.0) && (x < 2.0 * L / 3.0)) ? 1.0 : 0.0); }

static inline double
bool2 (const double x, const double L)
{ return (x <= L / 3.0 ? 1.0 : 0.0); }

static inline double
bool3 (const double x, const double L)
{ return (x >= 2.0 * L / 3.0 ? 1.0 : 0.0); }

static inline double
bool4 (const double x, const double L)
{ return (x < L / 4.0 ? 1.0 : 0.0); }

static inline double
bool5 (const double x, const double L)
{ return (x > 3.0 * L / 4.0 ? 1.0 : 0.0); }

static inline double
bool6 (const double x, const double L)
{ return (x >= L / 4.0 ? 1.0 : 0.0); }

static inline double
bool7 (const double x, const double L)
{ return (x <= 3.0 * L / 4.0 ? 1.0 : 0.0); }

double
doping  (double x, double y,
         double L, double H)
{
  constexpr double N_plus  = 1.0e25;
  constexpr double N_minus = 1.0e24;
  constexpr double P_plus  = 1.0e25;
  constexpr double P_minus = 1.0e24;

  double Na  = P_minus;
  Na += P_plus  * gaussian (y, -H / 2.0, H / 10.0) *
    (bool1 (x, L) +
     gaussian (x, L / 3.0, L / 10.0) * bool2 (x, L) +
     gaussian (x, 2.0 * L / 3.0, L / 10.0) * bool3 (x, L));

  double Nd  = 0.0;
  Nd += N_plus  * gaussian (y, 0, H / 4.0) *
    (bool4 (x, L) + bool5 (x, L) +
     gaussian (x, L / 4.0, L / 10.0) * bool6 (x, L) +
     gaussian (x, 3.0 * L / 4.0, L / 10.0) * bool7 (x, L));

  return (Nd - Na);
}

double
signedlog (double x)
{ return (asinh (x / 2.0) / log (10.0)); }

