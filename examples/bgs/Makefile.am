BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMIOLIB    = @OCTAVE_FILE_IO_LDFLAGS@ @OCTAVE_LINK_FLAGS@ @OCTAVE_FILE_IO_LIB@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_FILE_IO_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

bin_PROGRAMS = bgs_test

bgs_test_SOURCES = bgs_test.cpp
bgs_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB) $(BIMLISLIB) 
