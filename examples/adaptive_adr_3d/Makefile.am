BIMMUMPSLIB  = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB    = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB  = ../../src/libbimp4est.la @P4EST_LIB@
BIMLINALGLIB = ../../src/libbimlinalg.la 
BIMIOLIB     = @OCTAVE_FILE_IO_LDFLAGS@ @OCTAVE_LINK_FLAGS@ @OCTAVE_FILE_IO_LIB@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_FILE_IO_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

CLEANFILES = 

bin_PROGRAMS = p4est_adr_test_1_metrics p4est_dr_test_2_metrics

p4est_adr_test_1_metrics_SOURCES =  p4est_adr_test_1_metrics.cpp 
p4est_adr_test_1_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_adr_test_1_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_adr_test_1_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) @P4EST_LIB@
CLEANFILES += p4est_adr_test_1_metrics*.octbin.gz p4est_adr_test_1_metrics*.pvtu p4est_adr_test_1_metrics*.visit p4est_adr_test_1_metrics*.vtu

p4est_dr_test_2_metrics_SOURCES =  p4est_dr_test_2_metrics.cpp
p4est_dr_test_2_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_2_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_2_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) $(BIMLINALGLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_2_metrics*.octbin.gz p4est_dr_test_2_metrics*.pvtu p4est_dr_test_2_metrics*.visit p4est_dr_test_2_metrics*.vtu
