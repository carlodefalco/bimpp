#ifndef MOSFET_DOPING
#  define MOSFET_DOPING

#include <tmesh.h>

int
doping_driven_refinement (std::function<double (tmesh::idx_t, tmesh::idx_t)> p);

#endif
