BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMIOLIB    = ../../src/libbimio.la @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@
BIMLIB           = ../../src/libbim.la

AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

bin_PROGRAMS = p4est_refine_segment_test p4est_ref_test p4est_ref_lnodes_test p4est_palla_2d \
               p4est_palla_2d_flow simple_newton lis_csr_example simple simple_system cahn_hilliard cahn_hilliard_system
CLEANFILES =

p4est_refine_segment_test_SOURCES =  p4est_refine_segment_test.cpp
p4est_refine_segment_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_refine_segment_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_refine_segment_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
CLEANFILES += p4est_refine_segment_test.octbin*.gz p4est_refine_segment_test*.pvtu p4est_refine_segment_test*.visit p4est_refine_segment_test_*.vtu

p4est_ref_test_SOURCES =  p4est_ref_test.cpp mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
p4est_ref_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_ref_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_ref_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
CLEANFILES += p4est_ref_test.octbin*.gz p4est_ref_test*.pvtu p4est_ref_test*.visit p4est_ref_test_*.vtu

p4est_ref_lnodes_test_SOURCES =  p4est_ref_lnodes_test.cpp mosfet_doping_2d.cpp mosfet_connectivity_2d.cpp
p4est_ref_lnodes_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_ref_lnodes_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_ref_lnodes_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
CLEANFILES += p4est_ref_lnodes_test.octbin*.gz p4est_ref_lnodes_test*.pvtu p4est_ref_lnodes_test*.visit p4est_ref_lnodes_test_*.vtu

p4est_palla_2d_SOURCES =  p4est_palla_2d.cpp
p4est_palla_2d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_palla_2d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_palla_2d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
CLEANFILES += palla_initial_mesh*.pvtu palla_initial_mesh*.visit palla_initial_mesh*.vtu

p4est_palla_2d_flow_SOURCES =  p4est_palla_2d_flow.cpp
p4est_palla_2d_flow_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_palla_2d_flow_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_palla_2d_flow_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += palla_initial_mesh*.pvtu palla_initial_mesh*.visit palla_initial_mesh*.vtu

simple_newton_SOURCES =  simple_newton.cpp
simple_newton_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
simple_newton_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
simple_newton_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += 

simple_system_SOURCES =  simple_system.cpp
simple_system_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
simple_system_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
simple_system_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += 

cahn_hilliard_SOURCES =  cahn_hilliard.cpp
cahn_hilliard_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
cahn_hilliard_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
cahn_hilliard_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += 

cahn_hilliard_system_SOURCES =  cahn_hilliard_system.cpp
cahn_hilliard_system_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
cahn_hilliard_system_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
cahn_hilliard_system_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += 

simple_SOURCES =  simple.cpp
simple_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ -DBIM_TIMING
simple_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
simple_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += 


lis_csr_example_SOURCES =  lis_csr_example.cpp lis_csr_example_helper.cpp
lis_csr_example_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
lis_csr_example_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
lis_csr_example_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMLISLIB) @P4EST_LIB@
CLEANFILES += 
