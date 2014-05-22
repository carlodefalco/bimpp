BIMLIB = ../src/libbim.la ../src/libbimmumps.la
AM_CPPFLAGS =  -I$(top_srcdir)/src  @OCTAVE_COMPILE_FLAGS@ @MUMPS_CPPFLAGS@ 
AM_LDFLAGS = @OCTAVE_LINK_FLAGS@ @MUMPS_LDFLAGS@ 
AUTOMAKE_OPTIONS = nostdinc

check_PROGRAMS = sparse_test bim3_test

sparse_test_SOURCES = sparse_test.cpp
sparse_test_LDADD = $(BIMLIB)

bim3_test_SOURCES = bim3_test.cpp
bim3_test_LDADD = $(BIMLIB)

TESTS = sparse_test bim3_test
