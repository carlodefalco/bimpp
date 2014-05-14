AM_CPPFLAGS = -I$(top_srcdir)/src 
AUTOMAKE_OPTIONS = nostdinc

BIMXXHEADERS = mesh.h operators.h bim_sparse.h mumps_class.h mumps_solve.h 
BIMXXSRC = mesh.cpp operators.cpp bim_sparse.cpp
BIMXXMUMPSSRC = mumps_class.cpp mumps_solve.cpp
include_HEADERS = $(BIMXXHEADERS)

lib_LTLIBRARIES = libbim.la libbimmumps.la
libbim_la_SOURCES = $(BIMXXSRC)
libbim_la_LDFLAGS =  
libbim_la_CPPFLAGS = $(AM_CPPFLAGS)

libbimmumps_la_SOURCES = $(BIMXXMUMPSSRC)
libbimmumps_la_LDFLAGS = @MUMPS_LDFLAGS@
libbimmumps_la_CPPFLAGS = $(AM_CPPFLAGS) @MUMPS_CPPFLAGS@
libbimmumps_la_LIBADD = @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@









