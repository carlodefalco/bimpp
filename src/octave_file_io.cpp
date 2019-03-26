/*
  Copyright (C) 2011,2012,2012,2013,2014,2015,2016,2017 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/

#include <octave_file_io.h>

/*! \file octave_file_io.cpp
  \brief Classes and methods for files in Octave's (compressed) native binary file format
*/

// FIXME: This is a workaround for a problem introduced in Octave 4.2
#define HAVE_ZLIB

#include <fstream>
#include <iostream>
#include <octave/octave-config.h>
#include <octave/zfstream.h>

#include <octave/parse.h>
#include <octave/interpreter.h>

#include <octave/load-save.h>
#include <octave/ls-oct-binary.h>
#include <octave/oct-map.h>
#include <cstring>

#ifdef HAVE_OCTAVE_44
#  define OCT_ISEMPTY isempty
#else
#  define OCT_ISEMPTY is_empty
#endif

using namespace octave;

//---------------------------------------------------------------------
//                Singleton class
//---------------------------------------------------------------------

/// Singleton class providing an interface to Octave file I/O.
class octave_file_io_intf
{

public:

  octave_file_io_intf () 
    : filename ("")
  {
#ifdef HAVE_OCTAVE_44
    try
      {
        int status = interp.execute ();
        if (status != 0)
          {
            std::cerr << "creating embedded Octave interpreter failed!"
                      << std::endl;
            exit (status);
          }
      }
    catch (const octave::exit_exception& ex)
      {
        std::cerr << "Octave interpreter exited with status = "
                  << ex.exit_status () << std::endl;
        exit (ex.exit_status ());
      }
    catch (const octave::execution_exception&)
      {
        std::cerr << "error encountered in Octave evaluator!" << std::endl;
      }
#endif
  };
  
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

#ifdef HAVE_OCTAVE_44
  interpreter interp;
#endif
  std::fstream file;
  gzifstream gzifile;
  gzofstream gzofile;
  std::string filename;
  int current_mode;

};


//---------------------------------------------------------------------
//                Methods in the singleton class
//---------------------------------------------------------------------
const load_save_format format = LS_BINARY;
static octave::mach_info::float_format flt_fmt = octave::mach_info::flt_fmt_unknown;
static bool swap = false;

static bool
check_gzip_magic (const std::string &);

static bool 
fexists (const std::string &);

int
octave_file_io_intf::fopen (const char *fname, std::ios::openmode m)
{

  if (file.is_open ())
    file.close ();
  
  filename = fname;
  file.open (fname, m);
  
  if (m & std::ios::in)
    {
      current_mode = 0;
      if (read_binary_file_header (file, swap, flt_fmt, true) != 0)
	return -1;
    }
  else if (m & (std::ios::out & ! std::ios::app))
    {
      write_header (file, format);
      current_mode = 2;
    }
  else if (m & std::ios::app)
    {
      //write_header (file, format);
      current_mode = 2;
    }

  return 0;
}

int 
octave_file_io_intf::fclose (void)
{

  if (file.is_open ())
    file.close ();

  filename = "";
  return 0;
}

int
octave_file_io_intf::gzfopen (const char *fname, std::ios::openmode m)
{

  if (file.is_open ())
    file.close ();
   
  if (m & std::ios::in)
    {
      gzifile.open (fname, m);
      filename = fname;
      current_mode = 1;
      if (read_binary_file_header (gzifile, swap, flt_fmt, true) != 0)
	return -1;
    }
  else if (m & std::ios::out)
    {
      gzofile.open (fname, m);
      filename = fname;
      write_header (gzofile, format);
      current_mode = 3;
    }
  else if (m & std::ios::app)
    {
      gzofile.open (fname, m);
      filename = fname;
      //write_header (gzofile, format);
      current_mode = 3;
    }
  else
    return -1;

  return 0;
}

int 
octave_file_io_intf::gzfclose (void)
{

  if (gzifile.is_open ())
    {
      gzifile.close ();
      filename = "";
    }

  if (gzofile.is_open ())
    {
      gzofile.close ();
      filename = "";
    }

  return 0;
}


int
octave_file_io_intf::read
(const std::string &varname)
{
  string_vector argv (1);
#ifndef HAVE_OCTAVE_44
  install_types ();
#endif
  argv(0) = varname;
  
  file.clear ();
  file.seekg (0);
  if (read_binary_file_header (file, swap, flt_fmt, true) != 0)
    return -1;

  octave_scalar_map m = do_load (file, filename.c_str (), format, 
				 flt_fmt, false, swap, true, argv, 
				 0, 1, 1).scalar_map_value ();
  
  buffer = m.contents (varname);
  if (error_state)
    return -1;
  else
    return 0;
}

int
octave_file_io_intf::gzread 
(const std::string &varname)
{

  string_vector argv (1);
#ifndef HAVE_OCTAVE_44
  install_types ();
#endif
  argv(0) = varname;
  
  gzifile.clear ();
  gzifile.seekg (0);
  if (read_binary_file_header (gzifile, swap, flt_fmt, true) != 0)
    return -1;

  octave_value v = do_load (gzifile, filename.c_str (), format, 
                            flt_fmt, false, swap, true, argv,
                            0, 1, 1);

  if (! v.OCT_ISEMPTY ())
    buffer =  v.scalar_map_value ().contents (varname);
  
  if (v.OCT_ISEMPTY () || error_state)
    return -1;
  else
    return 0;  
}

int
octave_file_io_intf::do_read 
(const std::string &varname)
{
  switch (current_mode)
    {
    case 2:
    case 3:
      return -1;
    case 0:
      return read (varname);
      break;
    case 1:
      return gzread (varname);
      break;
    default:
      return -1;
    }
}

int
octave_file_io_intf::do_write
(const std::string &varname)
{
  switch (current_mode)
    {
    case 0:
    case 1:
      return -1;
    case 2:
      return write (varname);
      break;
    case 3:
      return gzwrite (varname);
      break;
    default:
      return -1;
    }
}

int
octave_file_io_intf::write 
(const std::string &varname)
{
  if (! save_binary_data (file, buffer, varname, "", false, false))
    return -1;
  return 0;
}

int
octave_file_io_intf::gzwrite
(const std::string &varname)
{
  if (! save_binary_data (gzofile, buffer, varname, "", false, false))
    return -1;
  return 0;  
}

// void 
// octave_file_io_intf::get_data_shape (int* ishape) 
// {
//   dim_vector d = buffer.dims ();
//   for (int ii = 0; ii < get_data_rank (); ++ii)
//     *(ishape + ii) = d(ii);
// }
// void 
// octave_file_io_intf::set_data_shape (const int rank, const int* ishape) 
// {
//   dim_vector d;
//   if (rank > 1)
//     d.resize (rank);
//   else
//     {d.resize (2); d(0) = 1; d(1) = 1;}
//   for (int ii = 0; ii < rank; ++ii)
//     d(ii) = *(ishape + ii);	
//   buffer.resize (d);
// }



//---------------------------------------------------------------------
//       Main API Functions 
//---------------------------------------------------------------------

static octave_file_io_intf io_intf;

int
octave_io_open (const char* fname, const octave_io_mode mode_in, 
                octave_io_mode *mode_out)
{

  *mode_out = mode_in;
  
  std::ios::openmode imode = std::ios::in | std::ios::binary;
  std::ios::openmode omode = std::ios::binary;
  
  bool gziped = false;
  if (fexists (std::string (fname)))
    {
      gziped = check_gzip_magic (fname);
      omode |= std::ios::ate | std::ios::app;
    }
  else
    omode |= std::ios::out;

  int retval = 0;
  if (gziped) 
    {
      if (*mode_out == read_mode)
         (*mode_out = gz_read_mode);
      else
        if (*mode_out == write_mode)
          (*mode_out = gz_write_mode);
    }

  //std::cout << "mode_out = " << *mode_out << std::endl;
  switch (*mode_out)
    {
    case 0: // read uncompressed
      retval = io_intf.fopen (fname, imode);
      break;
    case 1: // read compressed
      retval = io_intf.gzfopen (fname, imode);
      break;
    case 2: // write uncompressed
      retval = io_intf.fopen (fname, omode);
      break;
    case 3: // write compressed
      retval = io_intf.gzfopen (fname, omode);
      break;
    default:
      retval = -1;
    }

  return (retval);
}

int
octave_io_close (void)
{
  io_intf.fclose ();
  io_intf.gzfclose ();
  return 0;
}

int
octave_load (const char* varname, octave_value &data)
{
  int res = 0;

  res = io_intf.do_read (varname);
  data = io_intf.buffer;

  return res;
}

int
octave_save (const char* varname, const octave_value &data)
{
  int res = 0;

  io_intf.set_data (data);

  res = io_intf.do_write (varname);
  
  io_intf.clear_data ();

  return res;
}

int
octave_clear (void)
{
  int res = 0;
  io_intf.clear_data ();
  return res;
}


//---------------------------------------------------------------------
//                    Static utility functions
//---------------------------------------------------------------------

static bool
check_gzip_magic (const std::string& fname)
{
  bool retval = false;
  std::ifstream file (fname.c_str ());
  OCTAVE_LOCAL_BUFFER (unsigned char, magic, 2);

  if (file.read (reinterpret_cast<char *> (magic), 2) 
      && magic[0] == 0x1f 
      && magic[1] == 0x8b)

    retval = true;

  file.close ();
  return retval;
}


static bool 
fexists (const std::string& fname)
{
  std::ifstream ifile (fname.c_str ());
  return (ifile.good ());
}
