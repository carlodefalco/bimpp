BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMIOLIB    = ../../src/libbimio.la @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

bin_PROGRAMS = p4est_refine_segment_test p4est_ref_test


p4est_refine_segment_test_SOURCES =  p4est_refine_segment_test.cpp
p4est_refine_segment_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
p4est_refine_segment_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_refine_segment_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@

p4est_ref_test_SOURCES =  p4est_ref_test.cpp mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
p4est_ref_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
p4est_ref_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_ref_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
