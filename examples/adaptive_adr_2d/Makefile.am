BIMMUMPSLIB = ../../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB   = ../../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMP4ESTLIB = ../../src/libbimp4est.la @P4EST_LIB@
BIMIOLIB    = ../../src/libbimio.la @OCTAVE_LINK_FLAGS@ @OCTAVE_LINK_ADD@

BIMLIB           = ../../src/libbim.la
AM_CPPFLAGS      =  -I$(top_srcdir)/src -I$(top_srcdir)/test -I$(top_builddir)/src @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_LINK_FLAGS@ 
AM_LDFLAGS       = @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

CLEANFILES = 

bin_PROGRAMS = p4est_operator_segment_test \
               p4est_operator_test \
               p4est_operator_test_2 p4est_operator_test_3 \
               p4est_estimator_test_1 p4est_estimator_test_1_uniform \
               p4est_adr_test_1_marker p4est_adr_test_1_metrics p4est_adr_test_1_uniform \
               p4est_adr_test_2_marker p4est_adr_test_2_metrics p4est_adr_test_2_isogeo \
               p4est_metrics_test_1 \
               p4est_dr_test_1_marker p4est_dr_test_1_metrics p4est_dr_test_1_uniform \
               p4est_dr_test_2_marker p4est_dr_test_2_metrics \
               p4est_dr_test_3_marker p4est_dr_test_3_metrics p4est_dr_test_3_uniform \
               p4est_dr_test_4_marker p4est_dr_test_4_metrics p4est_dr_test_4_uniform \
               p4est_dr_test_5_marker p4est_dr_test_5_metrics p4est_dr_test_5_uniform \
               p4est_test_discontinuous

p4est_operator_test_SOURCES =  p4est_operator_test.cpp
p4est_operator_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_operator_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_operator_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_operator_test*.octbin.gz p4est_operator_test.pvtu p4est_operator_test.visit p4est_operator_test_*.vtu

p4est_operator_segment_test_SOURCES =  p4est_operator_segment_test.cpp
p4est_operator_segment_test_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_operator_segment_test_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_operator_segment_test_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_operator_segment_test*.octbin.gz p4est_operator_segment_test.pvtu p4est_operator_segment_test.visit p4est_operator_segment_test*.vtu

p4est_operator_test_3_SOURCES =  p4est_operator_test_3.cpp 
p4est_operator_test_3_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_operator_test_3_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_operator_test_3_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_operator_test_3*.octbin.gz p4est_operator_test_3.pvtu p4est_operator_test_3.visit p4est_operator_test_3_*.vtu

p4est_operator_test_2_SOURCES =  p4est_operator_test_2.cpp
p4est_operator_test_2_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_operator_test_2_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_operator_test_2_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_operator_test_2*.octbin.gz p4est_operator_test_2.pvtu p4est_operator_test_2.visit p4est_operator_test_2_*.vtu

p4est_estimator_test_1_SOURCES =  p4est_estimator_test_1.cpp 
p4est_estimator_test_1_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_estimator_test_1_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_estimator_test_1_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_estimator_test_1*.octbin.gz p4est_estimator_test_1*.pvtu p4est_estimator_test_1*.visit p4est_estimator_test_1*.vtu

p4est_estimator_test_1_uniform_SOURCES =  p4est_estimator_test_1_uniform.cpp 
p4est_estimator_test_1_uniform_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_estimator_test_1_uniform_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_estimator_test_1_uniform_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_estimator_test_1*.octbin.gz p4est_estimator_test_1*.pvtu p4est_estimator_test_1*.visit p4est_estimator_test_1*.vtu

p4est_adr_test_1_marker_SOURCES =  p4est_adr_test_1_marker.cpp
p4est_adr_test_1_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_adr_test_1_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_adr_test_1_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_adr_test_1_marker*.octbin.gz p4est_adr_test_1_marker*.pvtu p4est_adr_test_1_marker*.visit p4est_adr_test_1_marker*.vtu

p4est_adr_test_1_metrics_SOURCES =  p4est_adr_test_1_metrics.cpp
p4est_adr_test_1_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_adr_test_1_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_adr_test_1_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_adr_test_1_metrics*.octbin.gz p4est_adr_test_1_metrics*.pvtu p4est_adr_test_1_metrics*.visit p4est_adr_test_1_metrics*.vtu

p4est_adr_test_1_uniform_SOURCES =  p4est_adr_test_1_uniform.cpp
p4est_adr_test_1_uniform_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_adr_test_1_uniform_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_adr_test_1_uniform_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_adr_test_1_uniform*.octbin.gz p4est_adr_test_1_uniform*.pvtu p4est_adr_test_1_uniform*.visit p4est_adr_test_1_uniform*.vtu

p4est_adr_test_2_marker_SOURCES =  p4est_adr_test_2_marker.cpp 
p4est_adr_test_2_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_adr_test_2_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_adr_test_2_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_adr_test_2_marker*.octbin.gz p4est_adr_test_2_marker*.pvtu p4est_adr_test_2_marker*.visit p4est_adr_test_2_marker*.vtu

p4est_adr_test_2_metrics_SOURCES =  p4est_adr_test_2_metrics.cpp 
p4est_adr_test_2_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_adr_test_2_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_adr_test_2_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_adr_test_2_metrics*.octbin.gz p4est_adr_test_2_metrics*.pvtu p4est_adr_test_2_metrics*.visit p4est_adr_test_2_metrics*.vtu

p4est_adr_test_2_isogeo_SOURCES =  p4est_adr_test_2_isogeo.cpp 
p4est_adr_test_2_isogeo_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_adr_test_2_isogeo_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_adr_test_2_isogeo_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_adr_test_2_isogeo*.octbin.gz p4est_adr_test_2_isogeo*.pvtu p4est_adr_test_2_isogeo*.visit p4est_adr_test_2_isogeo*.vtu

p4est_metrics_test_1_SOURCES =  p4est_metrics_test_1.cpp 
p4est_metrics_test_1_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_metrics_test_1_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_metrics_test_1_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_metrics_test_1*.octbin.gz p4est_metrics_test_1*.pvtu p4est_metrics_test_1*.visit p4est_metrics_test_1*.vtu

p4est_dr_test_1_marker_SOURCES =  p4est_dr_test_1_marker.cpp 
p4est_dr_test_1_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_1_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_1_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_1_marker*.octbin.gz p4est_dr_test_1_marker*.pvtu p4est_dr_test_1_marker*.visit p4est_dr_test_1_marker*.vtu

p4est_dr_test_1_metrics_SOURCES =  p4est_dr_test_1_metrics.cpp 
p4est_dr_test_1_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_1_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_1_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_1_metrics*.octbin.gz p4est_dr_test_1_metrics*.pvtu p4est_dr_test_1_metrics*.visit p4est_dr_test_1_metrics*.vtu

p4est_dr_test_1_uniform_SOURCES =  p4est_dr_test_1_uniform.cpp 
p4est_dr_test_1_uniform_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_1_uniform_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_1_uniform_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_1_uniform*.octbin.gz p4est_dr_test_1_uniform*.pvtu p4est_dr_test_1_uniform*.visit p4est_dr_test_1_uniform*.vtu

p4est_dr_test_2_marker_SOURCES =  p4est_dr_test_2_marker.cpp 
p4est_dr_test_2_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_2_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_2_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_2_marker*.octbin.gz p4est_dr_test_2_marker*.pvtu p4est_dr_test_2_marker*.visit p4est_dr_test_2_marker*.vtu

p4est_dr_test_2_metrics_SOURCES =  p4est_dr_test_2_metrics.cpp 
p4est_dr_test_2_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_2_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_2_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_2_metrics*.octbin.gz p4est_dr_test_2_metrics*.pvtu p4est_dr_test_2_metrics*.visit p4est_dr_test_2_metrics*.vtu

p4est_dr_test_3_marker_SOURCES =  p4est_dr_test_3_marker.cpp 
p4est_dr_test_3_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_3_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_3_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_3_marker*.octbin.gz p4est_dr_test_3_marker*.pvtu p4est_dr_test_3_marker*.visit p4est_dr_test_3_marker*.vtu

p4est_dr_test_3_metrics_SOURCES =  p4est_dr_test_3_metrics.cpp 
p4est_dr_test_3_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_3_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_3_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_3_metrics*.octbin.gz p4est_dr_test_3_metrics*.pvtu p4est_dr_test_3_metrics*.visit p4est_dr_test_3_metrics*.vtu

p4est_dr_test_3_uniform_SOURCES =  p4est_dr_test_3_uniform.cpp 
p4est_dr_test_3_uniform_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_3_uniform_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_3_uniform_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_3_uniform*.octbin.gz p4est_dr_test_3_uniform*.pvtu p4est_dr_test_3_uniform*.visit p4est_dr_test_3_uniform*.vtu

p4est_dr_test_4_marker_SOURCES =  p4est_dr_test_4_marker.cpp 
p4est_dr_test_4_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_4_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_4_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_4_marker*.octbin.gz p4est_dr_test_4_marker*.pvtu p4est_dr_test_4_marker*.visit p4est_dr_test_4_marker*.vtu

p4est_dr_test_4_metrics_SOURCES =  p4est_dr_test_4_metrics.cpp 
p4est_dr_test_4_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_4_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_4_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_4_metrics*.octbin.gz p4est_dr_test_4_metrics*.pvtu p4est_dr_test_4_metrics*.visit p4est_dr_test_4_metrics*.vtu

p4est_dr_test_4_uniform_SOURCES =  p4est_dr_test_4_uniform.cpp 
p4est_dr_test_4_uniform_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_4_uniform_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_4_uniform_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_4_uniform*.octbin.gz p4est_dr_test_4_uniform*.pvtu p4est_dr_test_4_uniform*.visit p4est_dr_test_4_uniform*.vtu

p4est_dr_test_5_marker_SOURCES =  p4est_dr_test_5_marker.cpp 
p4est_dr_test_5_marker_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_5_marker_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_5_marker_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_5_marker*.octbin.gz p4est_dr_test_5_marker*.pvtu p4est_dr_test_5_marker*.visit p4est_dr_test_5_marker*.vtu

p4est_dr_test_5_metrics_SOURCES =  p4est_dr_test_5_metrics.cpp 
p4est_dr_test_5_metrics_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_5_metrics_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_5_metrics_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_5_metrics*.octbin.gz p4est_dr_test_5_metrics*.pvtu p4est_dr_test_5_metrics*.visit p4est_dr_test_5_metrics*.vtu

p4est_dr_test_5_uniform_SOURCES =  p4est_dr_test_5_uniform.cpp 
p4est_dr_test_5_uniform_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_dr_test_5_uniform_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_dr_test_5_uniform_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_dr_test_5_uniform*.octbin.gz p4est_dr_test_5_uniform*.pvtu p4est_dr_test_5_uniform*.visit p4est_dr_test_5_uniform*.vtu

p4est_test_discontinuous_SOURCES =  p4est_test_discontinuous.cpp
p4est_test_discontinuous_CPPFLAGS =  $(AM_CPPFLAGS) @P4EST_CPPFLAGS@
p4est_test_discontinuous_LDFLAGS =  $(AM_LDFLAGS)  @P4EST_LDFLAGS@
p4est_test_discontinuous_LDADD = $(BIMLIB) $(BIMIOLIB) $(BIMP4ESTLIB) $(BIMMUMPSLIB) @P4EST_LIB@
CLEANFILES += p4est_test_discontinuous*.octbin.gz p4est_test_discontinuous*.pvtu p4est_test_discontinuous*.visit p4est_test_discontinuous*.vtu
