
#ifndef BIM_ORDERING_H
#  define BIM_ORDERING_H

#include <functional>

/// Ordering.
using ordering = std::function<size_t (p4est_gloidx_t)>;

template<size_t ntot = 1, size_t n = 0>
size_t
dof_ordering (p4est_gloidx_t gt)
{ return ntot*gt+n; }

template<>
size_t
dof_ordering<1,0> (p4est_gloidx_t gt);

extern ordering default_ord;
extern ordering ord20;
extern ordering ord21;
extern ordering ord30;
extern ordering ord31;
extern ordering ord32;

#endif
