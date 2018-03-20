#include <p8est.h>

constexpr p4est_topidx_t simple_conn_num_vertices = 8;
constexpr p4est_topidx_t simple_conn_num_trees = 1;
const double simple_conn_p[simple_conn_num_vertices*3] = 
  {0., 0., 0., 1., 0., 0., 0., 1., 0., 1., 1., 0.,
   0., 0., 1., 1., 0., 1., 0., 1., 1., 1., 1., 1.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*9] = 
  {1, 2, 3, 4, 5, 6, 7, 8, 1};
