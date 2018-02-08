/*
  Copyright (C) 2011,2012,2012,2013,2014,2015,2016,2017 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#ifndef OCTAVE_FILE_IO_H
# define OCTAVE_FILE_IO_H

#include <octave/oct.h>
#include <octave/octave.h>

//---------------------------------------------------------------------
//                Singleton class
//---------------------------------------------------------------------

// Do not export the singleton class definition
// and other symbols that are not to be used elsewhere.

//---------------------------------------------------------------------
//                API Functions
//---------------------------------------------------------------------

enum octave_io_mode 
  {
    read_mode = 0,
    gz_read_mode = 1,
    write_mode = 2,
    gz_write_mode = 3
  };

int
octave_io_open (const char*, const octave_io_mode, octave_io_mode*);

int
octave_io_close (void);

int 
octave_load (const char*, octave_value&);

int 
octave_save (const char*, const octave_value&);

int
octave_clear (void);


#endif
