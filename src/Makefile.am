AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src 
AUTOMAKE_OPTIONS = nostdinc

BIMXXHEADERS =  bimutils.h bim_timing.h mesh.h operators.h \
bim_sparse.h bim_sparse_distributed.h linear_solver.h mumps_class.h mumps_solve.h \
lis_class.h nonlinear_solver.h adaptive_inexact_newton_class.h \
backtracking_inexact_newton_class.h abstract_nonlinear_problem.h \
abstract_forcing_term.h forcing_class.h bgs_class.h schur_class.h rre.h \
tmesh.h bim_ordering.h quad_operators.h tmesh_3d.h quad_operators_3d.h \
bim_distributed_vector.h mpi_alltoallv2.h system_solver.h lis_class_distributed.h

BIMXXSRC = mesh.cpp operators.cpp adaptive_inexact_newton_class.cpp \
backtracking_inexact_newton_class.cpp forcing_class.cpp  #schur_class.cpp

BIMXXLINALGSRC = bim_sparse.cpp bim_sparse_distributed.cpp bim_distributed_vector.cpp \
mpi_alltoallv2.cpp bimutils.cpp bgs_class.cpp rre.cpp
BIMXXMUMPSSRC = mumps_class.cpp mumps_solve.cpp
BIMXXP4ESTSRC = tmesh.cpp bim_ordering.cpp quad_operators.cpp tmesh_3d.cpp quad_operators_3d.cpp
BIMXXLISSRC = lis_class.cpp lis_class_master.cpp lis_class_slaves.cpp lis_class_distributed.cpp
BIMXXMUMPSP4ESTSRC = system_solver.cpp
include_HEADERS = $(BIMXXHEADERS)

lib_LTLIBRARIES = libbim.la libbimlinalg.la libbimp4est.la libbimmumps.la libbimlis.la libbimsystemsolver.la
libbim_la_SOURCES = $(BIMXXSRC)
libbim_la_LDFLAGS = @OCTAVE_LINK_FLAGS@ @OCTAVE_FILE_IO_LDFLAGS@ @OCTAVE_FILE_IO_LIB@
libbim_la_CPPFLAGS = @OCTAVE_COMPILE_FLAGS@ @OCTAVE_FILE_IO_CPPFLAGS@ $(AM_CPPFLAGS)

libbimlinalg_la_SOURCES = $(BIMXXLINALGSRC)
libbimlinalg_la_CPPFLAGS = $(AM_CPPFLAGS) 
libbimlinalg_la_LIBADD = 

libbimmumps_la_SOURCES = $(BIMXXMUMPSSRC)
libbimmumps_la_LDFLAGS = #@MUMPS_LDFLAGS@ 
libbimmumps_la_CPPFLAGS = $(AM_CPPFLAGS) @MUMPS_CPPFLAGS@
libbimmumps_la_LIBADD =  libbimlinalg.la 

libbimlis_la_SOURCES = $(BIMXXLISSRC)
libbimlis_la_LDFLAGS = #@LIS_LDFLAGS@ 
libbimlis_la_CPPFLAGS = $(AM_CPPFLAGS) @LIS_CPPFLAGS@
libbimlis_la_LIBADD = libbimlinalg.la #@LIS_LIB@

libbimp4est_la_SOURCES = $(BIMXXP4ESTSRC)
libbimp4est_la_LDFLAGS = @P4EST_LDFLAGS@  @OCTAVE_LINK_FLAGS@
libbimp4est_la_CPPFLAGS = $(AM_CPPFLAGS) @P4EST_CPPFLAGS@ @OCTAVE_COMPILE_FLAGS@ @OCTAVE_FILE_IO_CPPFLAGS@
libbimp4est_la_LIBADD = libbimlinalg.la @P4EST_LIB@ 

libbimsystemsolver_la_SOURCES = $(BIMXXMUMPSP4ESTSRC)
libbimsystemsolver_la_LDFLAGS = @MUMPS_LDFLAGS@ @P4EST_LDFLAGS@
libbimsystemsolver_la_CPPFLAGS = $(AM_CPPFLAGS) @MUMPS_CPPFLAGS@ @P4EST_CPPFLAGS@
libbimsystemsolver_la_LIBADD = libbimlinalg.la @P4EST_LIB@ #@MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@ 
