#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([libbim], [0.0.1], [])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/mesh.cpp])
AC_CONFIG_HEADERS([src/config.h])

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

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/Makefile])
AC_OUTPUT
