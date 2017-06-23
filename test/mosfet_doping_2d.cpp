#include <cmath>
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

static inline double
doping  (const double x, const double y,
         const double L, const double H)
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

static inline double
signedlog (double x)
{ return (asinh (x / 2.0) / log (10.0)); }

int
doping_driven_refinement (std::function<double (tmesh:idx_t, tmesh:idx_t)> p);
{

  constexpr double L = 3.0e-6;
  constexpr double H = 1.0e-5;

  double maxy = 0, miny = 0, y = 0;
  double top = p(1, 0);

  maxy = miny = y = signedlog (doping (p(1, 0), p(1, 1), L, H));

  for (int ii = 1; ii < 4; ++ii)
    {
      y = signedlog (doping (p(0, ii), p(1, ii), L, H));
      maxy = maxy < y ? y : maxy;
      miny = miny > y ? y : miny;
      top  = top < p(1, ii) ? p(1, ii) : top;
    }

  double delta = maxy - miny;
  return ((top <= 0 && delta > .1) ? 1 : 0);
  
}
