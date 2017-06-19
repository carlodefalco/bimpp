/*! \file tmesh.h
  \brief Interface for p4est library
*/

/*
 * tmesh is the p4est interface for the user: it contains all the necessary 
 * methods to allocate a p4est, refine or coarsen it, and access the mesh 
 * through quadrant iteration.
 * It also has methods for input/output in various formats.
 */

#ifndef TMESH_H
#define TMESH_H


#include <octave_file_io.h>
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_vtk.h>
#include <p4est_mesh.h>
#include <p4est_extended.h>
#include <p4est_lnodes.h>

class
tmesh
{

public:

  using idx_t = unsigned int;
  
  class
  quadrant_t
  {
    double
    p (idx_t i, idx_t j);

    idx_t
    t (idx_t i);

    bool
    is_hanging (idx_t i);

    void
    get_parents (idx_t i, std::vector<idx_t> pv);

    idx_t
    e (idx_t i);
  };

  class
  quadrant_iterator : public *quadrant_t
  {
  public:
    void 
    operator++ ();
  };

  quadrant_iterator
  begin_quadrant_sweep ();
  
  quadrant_iterator
  return_quadrant_sweep ()
  { return static_cast<quadrant_iterator> nullptr; };
  
};


#endif /* TMESH_H */

