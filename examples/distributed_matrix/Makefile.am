BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMLINALGLIB = ../../src/libbimlinalg.la 
BIMIOLIB    = @OCTAVE_FILE_IO_LDFLAGS@ @OCTAVE_LINK_FLAGS@ @OCTAVE_FILE_IO_LIB@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_FILE_IO_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

bin_PROGRAMS = distributed_matrix_example distributed_vector_example

distributed_matrix_example_SOURCES =  distributed_matrix_example.cpp
distributed_matrix_example_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
distributed_matrix_example_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
distributed_matrix_example_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) @P4EST_LIB@

distributed_vector_example_SOURCES =  distributed_vector_example.cpp
distributed_vector_example_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
distributed_vector_example_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
distributed_vector_example_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMLINALGLIB) @P4EST_LIB@

