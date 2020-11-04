BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMLINALGLIB = ../../src/libbimlinalg.la 
BIMIOLIB    = @OCTAVE_FILE_IO_LDFLAGS@  @OCTAVE_FILE_IO_LIB@ @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_FILE_IO_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
BIMCHLIB         = ../../src/libbimsystemsolver.la
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

bin_PROGRAMS = swe_tg
CLEANFILES =

swe_tg_SOURCES =  swe_tg.cpp
swe_tg_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
swe_tg_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
swe_tg_LDADD = $(BIMLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) $(BIMIOLIB) @P4EST_LIB@
CLEANFILES +=
