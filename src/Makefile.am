AM_CPPFLAGS = -I$(top_srcdir)/src 
AUTOMAKE_OPTIONS = nostdinc

lib_LTLIBRARIES = libbim.la
include_HEADERS = mesh.h bim_sparse.h operators.h mumps_class.h mumps_solve.h
libbim_la_SOURCES = mesh.cpp bim_sparse.cpp operators.cpp mumps_class.cpp mumps_solve.cpp
libbim_la_LDFLAGS = @OCTAVE_LINK_FLAGS@ @MUMPS_LDFLAGS@
libbim_la_CPPFLAGS = $(AM_CPPFLAGS) @OCTAVE_COMPILE_FLAGS@ @MUMPS_CPPFLAGS@ 
libbim_la_LIBADD = @OCTAVE_LINK_ADD@ @MUMPS_LIBS@ @MUMPS_EXTRA_LIBS@ $(LDADD) 









