AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src 
AUTOMAKE_OPTIONS = nostdinc

BIMXXHEADERS = octave_file_io.h bim_timing.h mesh.h operators.h \
bim_sparse.h linear_solver.h mumps_class.h mumps_solve.h \
lis_class.h nonlinear_solver.h adaptive_inexact_newton_class.h \
backtracking_inexact_newton_class.h backtracking_inexact_newton_example8_class.h abstract_nonlinear_problem.h \
abstract_forcing_term.h forcing_class.h bgs_class.h schur_class.h rre.h \
tmesh.h quad_operators.h

BIMXXSRC = mesh.cpp operators.cpp bim_sparse.cpp adaptive_inexact_newton_class.cpp \
backtracking_inexact_newton_class.cpp backtracking_inexact_newton_example8_class.cpp forcing_class.cpp bgs_class.cpp rre.cpp \
 #schur_class.cpp

BIMXXIOSRC = octave_file_io.cpp
BIMXXMUMPSSRC = mumps_class.cpp mumps_solve.cpp
BIMXXP4ESTSRC = tmesh.cpp quad_operators.cpp
BIMXXLISSRC = lis_class.cpp lis_class_master.cpp lis_class_slaves.cpp
include_HEADERS = $(BIMXXHEADERS)

lib_LTLIBRARIES = libbim.la libbimio.la libbimmumps.la libbimlis.la libbimp4est.la  
libbim_la_SOURCES = $(BIMXXSRC)
libbim_la_LDFLAGS = @OCTAVE_LINK_FLAGS@  
libbim_la_CPPFLAGS = @OCTAVE_COMPILE_FLAGS@ $(AM_CPPFLAGS)

libbimio_la_SOURCES = $(BIMXXIOSRC)
libbimio_la_LDFLAGS = @OCTAVE_LINK_FLAGS@
libbimio_la_CPPFLAGS = @OCTAVE_COMPILE_FLAGS@ $(AM_CPPFLAGS)
libbimio_la_LIBADD = @OCTAVE_LINK_ADD@

libbimmumps_la_SOURCES = $(BIMXXMUMPSSRC)
libbimmumps_la_LDFLAGS = @MUMPS_LDFLAGS@
libbimmumps_la_CPPFLAGS = $(AM_CPPFLAGS) @MUMPS_CPPFLAGS@
libbimmumps_la_LIBADD = @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@

libbimlis_la_SOURCES = $(BIMXXLISSRC)
libbimlis_la_LDFLAGS = @LIS_LDFLAGS@
libbimlis_la_CPPFLAGS = $(AM_CPPFLAGS) @LIS_CPPFLAGS@
libbimlis_la_LIBADD = @LIS_LIB@

libbimp4est_la_SOURCES = $(BIMXXP4ESTSRC)
libbimp4est_la_LDFLAGS = @P4EST_LDFLAGS@ @OCTAVE_LINK_FLAGS@
libbimp4est_la_CPPFLAGS = $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@
libbimp4est_la_LIBADD = @P4EST_LIB@
