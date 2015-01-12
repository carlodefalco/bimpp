AM_CPPFLAGS = -I$(top_srcdir)/src 
AUTOMAKE_OPTIONS = nostdinc

BIMXXHEADERS = mesh.h operators.h bim_sparse.h linear_solver.h mumps_class.h mumps_solve.h \
lis_operators.h lis_class.h nonlinear_operators.h nonlinear_solvers.h nonlinear_examples.h
BIMXXSRC = mesh.cpp operators.cpp bim_sparse.cpp nonlinear_operators.cpp nonlinear_examples.cpp
BIMXXMUMPSSRC = mumps_class.cpp mumps_solve.cpp
BIMXXLISSRC = lis_operators.cpp lis_class.cpp
include_HEADERS = $(BIMXXHEADERS)

lib_LTLIBRARIES = libbim.la libbimmumps.la libbimlis.la
libbim_la_SOURCES = $(BIMXXSRC)
libbim_la_LDFLAGS =  
libbim_la_CPPFLAGS = $(AM_CPPFLAGS)

libbimmumps_la_SOURCES = $(BIMXXMUMPSSRC)
libbimmumps_la_LDFLAGS = @MUMPS_LDFLAGS@
libbimmumps_la_CPPFLAGS = $(AM_CPPFLAGS) @MUMPS_CPPFLAGS@
libbimmumps_la_LIBADD = @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@

libbimlis_la_SOURCES = $(BIMXXLISSRC)
libbimlis_la_LDFLAGS = @LIS_LDFLAGS@
libbimlis_la_CPPFLAGS = $(AM_CPPFLAGS) @LIS_CPPFLAGS@
libbimlis_la_LIBADD = @LIS_LIB@










