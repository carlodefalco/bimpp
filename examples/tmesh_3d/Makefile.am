BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMIOLIB    = ../../src/libbimio.la @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

bin_PROGRAMS = p4est_palla_3d p4est_palla_3d_flow
CLEANFILES = 

p4est_palla_3d_SOURCES =  p4est_palla_3d.cpp
p4est_palla_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_palla_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_palla_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
CLEANFILES += palla_initial_mesh*.pvtu palla_initial_mesh*.visit palla_initial_mesh*.vtu

p4est_palla_3d_flow_SOURCES =  p4est_palla_3d_flow.cpp
p4est_palla_3d_flow_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_palla_3d_flow_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_palla_3d_flow_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += palla_initial_mesh*.pvtu palla_initial_mesh*.visit palla_initial_mesh*.vtu
