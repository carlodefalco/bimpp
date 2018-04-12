/*
  Copyright (C) 2015,2016,2017 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <mpi.h>
#include <bim_timing.h>
#include <octave_file_io.h>
#include <cassert>

int main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  int rank, size;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  std::vector<octave_idx_type>      n;
  std::vector<double>   x;
  MPI_Status            status;

  if (rank == 0)
    {
      tic ();

      std::vector<MPI_Request> ireq (size - 1);
      std::vector<MPI_Request> xreq (size - 1);
      
      // resize vectors
      n.resize (10 * size);
      x.resize (10 * size);
      
      // receive data from other ranks
      // and fill non-local entries
      for (int irank = 1; irank < size; ++irank)
        {

          MPI_Irecv (&(n[0]) + irank * 10, 10, MPI_INT,
                     irank, 0, MPI_COMM_WORLD,
                     &(ireq[irank - 1]));      

          MPI_Irecv (&(x[0]) + irank * 10, 10, MPI_DOUBLE,
                     irank, 1, MPI_COMM_WORLD,
                     &(xreq[irank - 1]));  
        }
      
      // fill in local entries
      for (int ii = 0; ii < 10; ++ii)
        { n[ii] = ii; x[ii] = ii * 0.1413; };

      // check that all communications are completed
      for (int irank = 1; irank < size; ++irank)
        {
          MPI_Wait (&(ireq[irank - 1]), &status);
          MPI_Wait (&(xreq[irank - 1]), &status);
        }
      
      // save data to file
      ColumnVector oct_x (x.size (), 0.0);
      Array<octave_idx_type> oct_n (dim_vector (n.size (), 1), 0);

      std::copy_n (x.begin (), x.size (), oct_x.fortran_vec ());
      std::copy_n (n.begin (), n.size (), oct_n.fortran_vec ());

      octave_scalar_map the_map;
      the_map.assign ("x", oct_x);
      the_map.assign ("n", oct_n);
      
      octave_io_mode m = gz_write_mode;
      assert (octave_io_open ("tmp_io_test.octbin.gz", m, &m) == 0);
      assert (octave_save ("the_map", octave_value (the_map)) == 0);
      assert (octave_io_close () == 0);
      
      // load data from file and print
      octave_value tmp;
      m = gz_read_mode;
      assert (octave_io_open ("tmp_io_test.octbin.gz", m, &m) == 0);
      assert (octave_load ("the_map", tmp) == 0);
      ColumnVector x_new = tmp.scalar_map_value ().contents ("x").column_vector_value ();
      Array<octave_idx_type> n_new = tmp.scalar_map_value ().contents ("n").column_vector_value ();

      for (int ii = 0; ii < x_new.numel (); ++ii)
        assert (x_new(ii) == x[ii]);

      for (int ii = 0; ii < n_new.numel (); ++ii)
        assert (n_new(ii) == n[ii]);

      assert (octave_io_close () == 0);

        
      toc ("save data");
    }
  else
    {
      MPI_Request ireq;
      MPI_Request xreq;
      
      // resize vectors
      n.resize (10);
      x.resize (10);

      // fill in local entries
      for (int ii = 0; ii < 10; ++ii)
        { n[ii] = ii + 1000 * rank; x[ii] = ii * 0.1413 + 1000 * rank; };
      
      // send local entries to rank 0
      MPI_Isend (&(n[0]), 10, MPI_INT, 0,
                 0, MPI_COMM_WORLD, &ireq);

      MPI_Isend (&(x[0]), 10, MPI_DOUBLE, 0,
                 1, MPI_COMM_WORLD, &xreq);

      MPI_Wait (&(ireq), &status);
      MPI_Wait (&(xreq), &status);
    }

  MPI_Barrier (MPI_COMM_WORLD);


  MPI_Finalize ();
  return (0);
}
