#include <p8est.h>

constexpr p4est_topidx_t simple_conn_num_vertices = 12;
constexpr p4est_topidx_t simple_conn_num_trees = 2;
const double simple_conn_p[simple_conn_num_vertices*3] = 
  {0., 0., 0., 1., 0., 0., 0., 1., 0., 1., 1., 0.,
   0., 0., 1., 1., 0., 1., 0., 1., 1., 1., 1., 1.,
   2., 0., 0., 2., 1., 0., 2., 0., 1., 2., 1., 1.};
const p4est_topidx_t simple_conn_t[simple_conn_num_trees*9] = 
  {1, 2, 3, 4, 5, 6, 7, 8, 1,
   2, 9, 4, 10, 6, 11, 8, 12, 1};
