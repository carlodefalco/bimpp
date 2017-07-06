#ifndef HAVE_QUAD_OPERATORS_H
#define HAVE_QUAD_OPERATORS_H 1

#include "bim_sparse.h"
#include "operators.h"
#include "tmesh.h"


void bim2a_advection_diffusion (tmesh & mesh,
                                const std::vector<double> & alpha,
                                const std::vector<double> & psi,
                                sparse_matrix & A);

void bim2a_rhs (const tmesh & mesh,
                const std::vector<double> & f,
                const std::vector<double> & g,
                std::vector<double> & rhs);

#endif
