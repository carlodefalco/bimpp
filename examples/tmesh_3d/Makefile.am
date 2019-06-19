BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMIOLIB    = ../../src/libbimio.la @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

bin_PROGRAMS = p4est_palla_3d p4est_palla_3d_flow test_gradient_3d p4est_refine_box poisson_boltzmann
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

test_gradient_3d_SOURCES =  test_gradient_3d.cpp
test_gradient_3d_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
test_gradient_3d_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
test_gradient_3d_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
CLEANFILES += test_grad_3d_initial_mesh*.pvtu test_grad_3d_initial_mesh*.visit test_grad_3d_initial_mesh*.vtu

p4est_refine_box_SOURCES =  p4est_refine_box.cpp
p4est_refine_box_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_refine_box_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_refine_box_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_refine_box_*.pvtu p4est_refine_box_*.visit p4est_refine_box_*.vtu

poisson_boltzmann_SOURCES =  poisson_boltzmann.cpp
poisson_boltzmann_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
poisson_boltzmann_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
poisson_boltzmann_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) @P4EST_LIB@
CLEANFILES += poisson_boltzmann_initial_mesh*.pvtu poisson_boltzmann_initial_mesh*.visit poisson_boltzmann_mesh*.vtu
