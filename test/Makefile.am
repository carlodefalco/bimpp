BIMMUMPSLIB = ../src/libbimmumps.la @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@
BIMLISLIB = ../src/libbimlis.la @LIS_LIB@ @MUMPS_EXTRA_LIBS@
BIMLIB = ../src/libbim.la 
AM_CPPFLAGS =  -I$(top_srcdir)/src -I$(top_builddir)/src  @OCTAVE_COMPILE_FLAGS@ @MUMPS_CPPFLAGS@ @LIS_CPPFLAGS@
AM_LDFLAGS = @OCTAVE_LINK_FLAGS@ @MUMPS_LDFLAGS@ @LIS_LDFLAGS@
AUTOMAKE_OPTIONS = nostdinc

CLEANFILES = 

check_PROGRAMS = sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 bim3_stationary_test \
                 bim3_stationary_test2 bim3_stationary_test2_osc bim3_stationary_test3 \
                 bim3_stationary_test3_osc bim3_stationary_test4 bim3_timedependent_test \
                 bim3_timedependent_test2 bim3_timedependent_test2_osc bim3_timedependent_test3 \
                 bim3_timedependent_test3_osc 

if HAVE_LIS
  check_PROGRAMS += lis_test1 lis_stationary_test lis_timedependent_test \
                    bim3_nonlinear_test bim3_nonlinear_eq_test bim3_nonlinear_test2
endif
sparse_test_SOURCES = sparse_test.cpp
CLEANFILES += 
sparse_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_test_SOURCES = bim3_test.cpp
bim3_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_test2_SOURCES = bim3_test2.cpp
bim3_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_test3_SOURCES = bim3_test3.cpp
bim3_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_test4_SOURCES = bim3_test4.cpp
bim3_test4_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_stationary_test_SOURCES = bim3_stationary_test.cpp
bim3_stationary_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_stationary_test2_SOURCES = bim3_stationary_test2.cpp
bim3_stationary_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_stationary_test2_osc_SOURCES = bim3_stationary_test2_osc.cpp
bim3_stationary_test2_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_stationary_test3_SOURCES = bim3_stationary_test3.cpp
bim3_stationary_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_stationary_test3_osc_SOURCES = bim3_stationary_test3_osc.cpp
bim3_stationary_test3_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_stationary_test4_SOURCES = bim3_stationary_test4.cpp
bim3_stationary_test4_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_timedependent_test_SOURCES = bim3_timedependent_test.cpp
bim3_timedependent_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_timedependent_test2_SOURCES = bim3_timedependent_test2.cpp
bim3_timedependent_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_timedependent_test2_osc_SOURCES = bim3_timedependent_test2_osc.cpp
bim3_timedependent_test2_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_timedependent_test3_SOURCES = bim3_timedependent_test3.cpp
bim3_timedependent_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_timedependent_test3_osc_SOURCES = bim3_timedependent_test3_osc.cpp
bim3_timedependent_test3_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

if HAVE_LIS
  lis_test1_SOURCES = lis_test1.cpp
  lis_test1_LDADD = $(BIMLIB) $(BIMLISLIB)

  lis_stationary_test_SOURCES = lis_stationary_test.cpp
  lis_stationary_test_LDADD = $(BIMLIB) $(BIMLISLIB)

  lis_timedependent_test_SOURCES = lis_timedependent_test.cpp
  lis_timedependent_test_LDADD = $(BIMLIB) $(BIMLISLIB)

  bim3_nonlinear_eq_test_SOURCES = bim3_nonlinear_eq_test.cpp
  bim3_nonlinear_eq_test_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)

  bim3_nonlinear_test_SOURCES = bim3_nonlinear_test.cpp
  bim3_nonlinear_test_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)

  bim3_nonlinear_test2_SOURCES = bim3_nonlinear_test2.cpp
  bim3_nonlinear_test2_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)

endif

TESTS = sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 \
        bim3_stationary_test bim3_stationary_test2 bim3_stationary_test2_osc \
        bim3_stationary_test3 bim3_stationary_test3_osc bim3_stationary_test4 \
        bim3_timedependent_test bim3_timedependent_test2 bim3_timedependent_test2_osc \
        bim3_timedependent_test3 bim3_timedependent_test3_osc  

if HAVE_LIS
  TESTS+=lis_test1 lis_stationary_test lis_timedependent_test \
         bim3_nonlinear_eq_test bim3_nonlinear_test  bim3_nonlinear_test2
endif
