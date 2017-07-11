#ifndef HAVE_QUAD_OPERATORS_H
#define HAVE_QUAD_OPERATORS_H 1

#include "bim_sparse.h"
#include "operators.h"
#include "tmesh.h"

#include <functional>
#include <tuple>
#include <vector>

using func = std::function<double (const double &, const double &)>; // f(x, y).
using dirichlet_bcs = std::vector<std::tuple<const int &, const int &, const func &>>;

void bim2a_advection_diffusion (tmesh & mesh,
                                const std::vector<double> & alpha,
                                const std::vector<double> & psi,
                                sparse_matrix & A);

void bim2a_rhs (tmesh & mesh,
                const std::vector<double> & f,
                const std::vector<double> & g,
                std::vector<double> & rhs);

void bim2a_dirichlet_bc (tmesh & mesh, const dirichlet_bcs & bcs,
                         sparse_matrix & A, std::vector<double> & rhs);

#endif
