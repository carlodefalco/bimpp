/*
  Copyright (C) 2011,2012,2012,2013,2014,2015,2016,2017 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#ifndef OCTAVE_FILE_IO_H
# define OCTAVE_FILE_IO_H

// FIXME: This is a workaround for a problem introduced in Octave 4.2
#define HAVE_ZLIB

#include <fstream>
#include <octave/octave-config.h>
#include <octave/zfstream.h>

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>

#include <octave/load-save.h>
#include <octave/ls-oct-binary.h>
#include <octave/oct-map.h>
#include <cstring>

using namespace octave;

//---------------------------------------------------------------------
//                Singleton class
//---------------------------------------------------------------------

/// Singleton class providing an interface to Octave file I/O.
class octave_file_io_intf
{

public:

  octave_file_io_intf () 
    : filename ("") {};
  
  int fopen (const char *fname, std::ios::openmode m);
  int gzfopen (const char *fname, std::ios::openmode m);
  
  int fclose (void);
  int gzfclose (void);

  int do_read (const std::string &);
  int do_write (const std::string &);
  octave_value buffer;

  octave_value& get_data (void) { return buffer; };  
  void inline set_data (const octave_value& data) { buffer = data; };
  void inline clear_data (void) { buffer = 0; };

private:

  int read (const std::string &);
  int gzread (const std::string &);

  int write (const std::string &);
  int gzwrite (const std::string &);

  std::fstream file;
  gzifstream gzifile;
  gzofstream gzofile;
  std::string filename;
  int current_mode;

};

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
