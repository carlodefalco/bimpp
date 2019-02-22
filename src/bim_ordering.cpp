#include <p4est.h>
#include <bim_ordering.h>

template<>
size_t
dof_ordering<1,0> (p4est_gloidx_t gt)
{ return static_cast<size_t> (gt); }

ordering default_ord = dof_ordering<1,0>;
