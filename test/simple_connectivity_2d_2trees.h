#include <p4est.h>

constexpr p4est_topidx_t simple_conn_num_vertices = 6;
constexpr p4est_topidx_t simple_conn_num_trees = 2;
const double simple_conn_p[simple_conn_num_vertices*2] = {0., 0., 1., 0.,  1., 1., 0., 1., 2., 0., 2., 1.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*5] = {1, 2, 3, 4, 1, 2, 5, 6, 3, 1};
