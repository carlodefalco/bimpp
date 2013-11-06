# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([libbim], [0.0.1], [])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/mesh.cpp])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_MACRO_DIR([m4])

# Init libtool.
ifdef([LT_INIT], [], [
  errprint([error: you must have libtool 2.2.2 or a more recent version
])
  m4exit([1])])

LT_PREREQ([2.2.2])
LT_INIT([disable-static dlopen])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

AC_LANG_PUSH([C++])

# Checks for libraries.
AC_CHECK_LIB([mpi], [MPI_Init])
AC_CHECK_LIB([dmumps], [dmumps_c])

# Checks for header files.
AC_CHECK_HEADER(dmumps_c.h, ,AC_MSG_ERROR([dmumps_c.h not found.]))

AC_LANG_POP([C++])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T


# Checks for octave build script "mkoctfile".
AC_CHECK_PROG(HAVE_MKOCTFILE,mkoctfile,yes)
AM_CONDITIONAL(BUILD_OCT_TESTS,test x"${HAVE_MKOCTFILE}" == x"yes")
if [test x"${HAVE_MKOCTFILE}" == x"yes"]; then

# use mkoctfile to get info on octave flags
   AC_SUBST(OCTAVE_COMPILE_FLAGS, "$(env -i PATH=$PATH mkoctfile -p INCFLAGS) \
$(env -i PATH=$PATH mkoctfile -p CPPFLAGS) \
$(env -i PATH=$PATH mkoctfile -p CXXPICFLAG) \
$(env -i PATH=$PATH mkoctfile -p ALL_CXXFLAGS)")
   AC_SUBST(OCTAVE_LINK_FLAGS, "$(env -i PATH=$PATH mkoctfile -p CPPFLAGS) \
$(env -i PATH=$PATH mkoctfile -p ALL_CXXFLAGS) \
$(env -i PATH=$PATH mkoctfile -p RDYNAMIC_FLAG) \
$(env -i PATH=$PATH mkoctfile -p ALL_LDFLAGS)")
   AC_SUBST(OCTAVE_LINK_ADD, "$(env -i PATH=$PATH mkoctfile -p LFLAGS) \
$(env -i PATH=$PATH mkoctfile -p OCTAVE_LIBS) \
$(env -i PATH=$PATH mkoctfile -p LAPACK_LIBS) \
$(env -i PATH=$PATH mkoctfile -p BLAS_LIBS) \
$(env -i PATH=$PATH mkoctfile -p FFTW_LIBS) \
$(env -i PATH=$PATH mkoctfile -p LIBREADLINE) \
$(env -i PATH=$PATH mkoctfile -p LIBS) \
$(env -i PATH=$PATH mkoctfile -p FLIBS) \
")
else
   AC_MSG_WARN([You do not seem to have a complete install of Octave on your computer.])
fi

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 data/Makefile     
                 test/Makefile])
AC_OUTPUT
