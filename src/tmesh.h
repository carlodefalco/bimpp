/*
  Copyright (C) 2017 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

/*! \file tmesh.h
  \brief Interface for p4est library
*/

#ifndef TMESH_H
#define TMESH_H


#include <octave_file_io.h>

#include <mpi.h>

#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_vtk.h>

#include <cassert>
#include <functional>
#include <array>
#include <vector>


class
tmesh
{

public:

  using idx_t = unsigned int;
  
  class
  quadrant_t
  {

  public:

    quadrant_t (tmesh *_tmesh,
                p4est_topidx_t _tree = 0,
                p4est_quadrant_t *_quadrant = nullptr) :
      the_tmesh(_tmesh), the_tree(_tree), the_quadrant(_quadrant)
    {  };

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

    void
    update (p4est_topidx_t tree,
            p4est_quadrant_t *q);
    
  private:

    tmesh               *the_tmesh;
    p4est_topidx_t        the_tree;
    p4est_quadrant_t *the_quadrant;
    double                vxyz[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
  };

  class
  quadrant_iterator  
  {

  public:
    
    void 
    operator++ ();
    
    quadrant_t&
    operator* ()
    { return *(this->data); };
    
    const quadrant_t&
    operator* () const
    { return *(this->data); };

    const quadrant_t*
    get_data () const
    { return this->data; };

    quadrant_t*
    get_data ()
    { return this->data; };

    bool
    operator== (const quadrant_iterator& other)
    { return (this->get_data () == other.get_data ()); };
    
    bool
    operator!= (const quadrant_iterator& other)
    { return !((*this) == other); };
    
    quadrant_iterator (quadrant_t *_data = nullptr) :
      data (_data)
    { };
    
  private:
    quadrant_t *data;      
  };

  
  tmesh ()
    : p4est (nullptr), conn (nullptr), current_quadrant (this, 0, nullptr)
  {  };

  tmesh (const char *filename)
    : p4est (nullptr), current_quadrant (this, 0, nullptr)
  { read_connectivity (filename); };

  ~tmesh ();

  void
  read_connectivity (const char *filename,
                     int source = 0,
                     MPI_Comm comm = MPI_COMM_WORLD);

  quadrant_iterator
  begin_quadrant_sweep ();
  
  quadrant_iterator
  end_quadrant_sweep ()
  { return quadrant_iterator (); };

  void
  set_refine_marker
  (std::function<int (quadrant_iterator)> fun)
  { refine_marker = fun; };

  void
  set_coarsen_marker
  (std::function<int (quadrant_iterator[4])> fun)
  { coarsen_marker = fun; };

  void
  refine (int recursive = 0, int partforcoarsen = 0);

  // temporarily public untli the API is stable
  p4est_t              *p4est;
  p4est_connectivity_t *conn;
  quadrant_t            current_quadrant;
  
private:

  std::function<int (quadrant_iterator)> refine_marker;
  std::function<int (quadrant_iterator[4])> coarsen_marker;

  static int
  refine_callback (p4est_t*, p4est_topidx_t, p4est_quadrant_t*);

  static int
  coarsen_callback (p4est_t*, p4est_topidx_t, p4est_quadrant_t* []);

  void
  update_quadrant (p4est_topidx_t tree,
                   p4est_quadrant_t *q);
  
};


#endif /* TMESH_H */

