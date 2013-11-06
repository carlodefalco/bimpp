
include_HEADERS = mesh.h sparse.h operators.h mumps_class.h mumps_solve.h
lib_LTLIBRARIES = libbim.la
libbim_la_SOURCES = mesh.cpp sparse.cpp operators.cpp mumps_class.cpp mumps_solve.cpp
libbim_la_LIBADD = @OCTAVE_LINK_ADD@ mesh.o sparse.o mesh.o operators.o mumps_class.o mumps_solve.o libmesh.a

#libbim_a_LDFLAGS = @OCTAVE_LINK_FLAGS@ 
libbim_la_CPPFLAGS = @OCTAVE_COMPILE_FLAGS@ $(AM_CPPFLAGS)






