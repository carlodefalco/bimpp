BIMLIB = ../src/libbim.la ../src/libbimmumps.la
AM_CPPFLAGS =  -I$(top_srcdir)/src  @OCTAVE_COMPILE_FLAGS@ @MUMPS_CPPFLAGS@ 
AM_LDFLAGS = @OCTAVE_LINK_FLAGS@ @MUMPS_LDFLAGS@ 
AUTOMAKE_OPTIONS = nostdinc

check_PROGRAMS = sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 bim3_stationary_test bim3_stationary_test2 bim3_stationary_test2_osc bim3_stationary_test3 bim3_stationary_test3_osc bim3_stationary_test4 bim3_timedependent_test bim3_timedependent_test2 bim3_timedependent_test2_osc bim3_timedependent_test3 bim3_timedependent_test3_osc

sparse_test_SOURCES = sparse_test.cpp
sparse_test_LDADD = $(BIMLIB)

bim3_test_SOURCES = bim3_test.cpp
bim3_test_LDADD = $(BIMLIB)

bim3_test2_SOURCES = bim3_test2.cpp
bim3_test2_LDADD = $(BIMLIB)

bim3_test3_SOURCES = bim3_test3.cpp
bim3_test3_LDADD = $(BIMLIB)

bim3_test4_SOURCES = bim3_test4.cpp
bim3_test4_LDADD = $(BIMLIB)

bim3_stationary_test_SOURCES = bim3_stationary_test.cpp
bim3_stationary_test_LDADD = $(BIMLIB)

bim3_stationary_test2_SOURCES = bim3_stationary_test2.cpp
bim3_stationary_test2_LDADD = $(BIMLIB)

bim3_stationary_test2_osc_SOURCES = bim3_stationary_test2_osc.cpp
bim3_stationary_test2_osc_LDADD = $(BIMLIB)

bim3_stationary_test3_SOURCES = bim3_stationary_test3.cpp
bim3_stationary_test3_LDADD = $(BIMLIB)

bim3_stationary_test3_osc_SOURCES = bim3_stationary_test3_osc.cpp
bim3_stationary_test3_osc_LDADD = $(BIMLIB)

bim3_stationary_test4_SOURCES = bim3_stationary_test4.cpp
bim3_stationary_test4_LDADD = $(BIMLIB)

bim3_timedependent_test_SOURCES = bim3_timedependent_test.cpp
bim3_timedependent_test_LDADD = $(BIMLIB)

bim3_timedependent_test2_SOURCES = bim3_timedependent_test2.cpp
bim3_timedependent_test2_LDADD = $(BIMLIB)

bim3_timedependent_test2_osc_SOURCES = bim3_timedependent_test2_osc.cpp
bim3_timedependent_test2_osc_LDADD = $(BIMLIB)

bim3_timedependent_test3_SOURCES = bim3_timedependent_test3.cpp
bim3_timedependent_test3_LDADD = $(BIMLIB)

bim3_timedependent_test3_osc_SOURCES = bim3_timedependent_test3_osc.cpp
bim3_timedependent_test3_osc_LDADD = $(BIMLIB)

TESTS = sparse_test bim3_test bim3_test2 bim3_test3 bim3_test4 bim3_stationary_test bim3_stationary_test2 bim3_stationary_test2_osc bim3_stationary_test3 bim3_stationary_test3_osc bim3_stationary_test4 bim3_timedependent_test bim3_timedependent_test2 bim3_timedependent_test2_osc bim3_timedependent_test3 bim3_timedependent_test3_osc
