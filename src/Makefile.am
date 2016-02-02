AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src 
AUTOMAKE_OPTIONS = nostdinc

BIMXXHEADERS = bim_timing.h mesh.h operators.h bim_sparse.h linear_solver.h mumps_class.h mumps_solve.h \
lis_class.h nonlinear_solver.h adaptive_inexact_newton_class.h \
backtracking_inexact_newton_class.h abstract_nonlinear_problem.h \
abstract_forcing_term.h forcing_class.h bgs_class.h schur_class.h rre.h

BIMXXSRC = mesh.cpp operators.cpp bim_sparse.cpp adaptive_inexact_newton_class.cpp \
backtracking_inexact_newton_class.cpp forcing_class.cpp bgs_class.cpp rre.cpp \
schur_class.cpp

BIMXXMUMPSSRC = mumps_class.cpp mumps_solve.cpp
BIMXXLISSRC = lis_class.cpp lis_class_master.cpp lis_class_slaves.cpp
BIMXXBGSSRC = bgs_class.cpp rre.cpp
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
