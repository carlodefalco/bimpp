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
                 bim3_timedependent_test3_osc linear_solver_test

if HAVE_LIS
  check_PROGRAMS += lis_test1 lis_stationary_test lis_timedependent_test \
                    bim3_nonlinear_test bim3_nonlinear_eq_test bim3_nonlinear_test2 \
                    bim3_nonlinear_test4 bim3_nonlinear_test3
endif
sparse_test_SOURCES = sparse_test.cpp
CLEANFILES += first_solve.txt second_solve.txt
sparse_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)

bim3_test_SOURCES = bim3_test.cpp
bim3_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += SG.m mesh_out.m

bim3_test2_SOURCES = bim3_test2.cpp
bim3_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += SG2.m mesh_out.m

bim3_test3_SOURCES = bim3_test3.cpp
bim3_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += SG3.m mesh_out.m

bim3_test4_SOURCES = bim3_test4.cpp
bim3_test4_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += UP.m mesh_out.m

bim3_stationary_test_SOURCES = bim3_stationary_test.cpp
bim3_stationary_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution.txt

bim3_stationary_test2_SOURCES = bim3_stationary_test2.cpp
bim3_stationary_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution2.txt

bim3_stationary_test2_osc_SOURCES = bim3_stationary_test2_osc.cpp
bim3_stationary_test2_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution2_OSC.txt

bim3_stationary_test3_SOURCES = bim3_stationary_test3.cpp
bim3_stationary_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution3.txt

bim3_stationary_test3_osc_SOURCES = bim3_stationary_test3_osc.cpp
bim3_stationary_test3_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution3_OSC.txt

bim3_stationary_test4_SOURCES = bim3_stationary_test4.cpp
bim3_stationary_test4_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += solution4.txt

bim3_timedependent_test_SOURCES = bim3_timedependent_test.cpp
bim3_timedependent_test_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest1.txt

bim3_timedependent_test2_SOURCES = bim3_timedependent_test2.cpp
bim3_timedependent_test2_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest2.txt

bim3_timedependent_test2_osc_SOURCES = bim3_timedependent_test2_osc.cpp
bim3_timedependent_test2_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest2_osc.txt

bim3_timedependent_test3_SOURCES = bim3_timedependent_test3.cpp
bim3_timedependent_test3_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest2.txt

bim3_timedependent_test3_osc_SOURCES = bim3_timedependent_test3_osc.cpp
bim3_timedependent_test3_osc_LDADD = $(BIMLIB) $(BIMMUMPSLIB)
CLEANFILES += Solution_TimeTest3_osc.txt

if HAVE_LIS
  lis_test1_SOURCES = lis_test1.cpp
  lis_test1_LDADD = $(BIMLIB) $(BIMLISLIB)

  lis_stationary_test_SOURCES = lis_stationary_test.cpp
  lis_stationary_test_LDADD = $(BIMLIB) $(BIMLISLIB)
  CLEANFILES += LIS_Solution.txt

  lis_timedependent_test_SOURCES = lis_timedependent_test.cpp
  lis_timedependent_test_LDADD = $(BIMLIB) $(BIMLISLIB)
  CLEANFILES += Lis_Solution_TimeTest.txt

  bim3_nonlinear_eq_test_SOURCES = bim3_nonlinear_eq_test.cpp
  bim3_nonlinear_eq_test_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)
  CLEANFILES += Solution_Equation.txt

  bim3_nonlinear_test_SOURCES = bim3_nonlinear_test.cpp
  bim3_nonlinear_test_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)
  CLEANFILES += Lis_Solution_NonLinear.txt

  bim3_nonlinear_test2_SOURCES = bim3_nonlinear_test2.cpp
  bim3_nonlinear_test2_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)
  CLEANFILES += Lis_Solution_NonLinear2.txt

  bim3_nonlinear_test3_SOURCES = bim3_nonlinear_test3.cpp
  bim3_nonlinear_test3_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)
  CLEANFILES += Lis_Solution_NonLinear3.txt

  bim3_nonlinear_test4_SOURCES = bim3_nonlinear_test4.cpp
  bim3_nonlinear_test4_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)
  CLEANFILES += Lis_Solution_NonLinear4.txt

  linear_solver_test_SOURCES = linear_solver_test.cpp
  linear_solver_test_LDADD = $(BIMLIB) $(BIMLISLIB) $(BIMMUMPSLIB)
  CLEANFILES += MUMPS_solution.txt LIS_Solution.txt
endif

TESTS = sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 \
        bim3_stationary_test bim3_stationary_test2 bim3_stationary_test2_osc \
        bim3_stationary_test3 bim3_stationary_test3_osc bim3_stationary_test4 \
        bim3_timedependent_test bim3_timedependent_test2 bim3_timedependent_test2_osc \
        bim3_timedependent_test3 bim3_timedependent_test3_osc linear_solver_test

if HAVE_LIS
  TESTS+=lis_test1 lis_stationary_test lis_timedependent_test \
         bim3_nonlinear_eq_test bim3_nonlinear_test  bim3_nonlinear_test2 \
         bim3_nonlinear_test3 bim3_nonlinear_test4
endif
