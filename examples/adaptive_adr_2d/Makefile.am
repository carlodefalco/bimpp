BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMIOLIB    = ../../src/libbimio.la @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

bin_PROGRAMS = p4est_operator_segment_test \
	       p4est_operator_test_2 p4est_operator_test_3 

p4est_operator_segment_test_SOURCES =  p4est_operator_segment_test.cpp
p4est_operator_segment_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
p4est_operator_segment_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_operator_segment_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
##CLEANFILES += p4est_operator_segment_test*.octbin.gz p4est_operator_segment_test.pvtu p4est_operator_segment_test.visit p4est_operator_segment_test*.vtu

p4est_operator_test_3_SOURCES =  p4est_operator_test_3.cpp 
p4est_operator_test_3_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
p4est_operator_test_3_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_operator_test_3_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
##CLEANFILES += p4est_operator_test_3*.octbin.gz p4est_operator_test_3.pvtu p4est_operator_test_3.visit p4est_operator_test_3_*.vtu

p4est_operator_test_2_SOURCES =  p4est_operator_test_2.cpp
p4est_operator_test_2_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
p4est_operator_test_2_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_operator_test_2_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
##CLEANFILES += p4est_operator_test_2*.octbin.gz p4est_operator_test_2.pvtu p4est_operator_test_2.visit p4est_operator_test_2_*.vtu
