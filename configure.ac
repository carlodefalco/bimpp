# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_CONFIG_MACRO_DIR([m4])
AC_INIT([libbim], [0.0.1], [])
AC_CONFIG_SRCDIR([src/mesh.cpp])
AC_CONFIG_HEADERS([src/config.h])

AM_INIT_AUTOMAKE([subdir-objects])

# Init libtool.
LT_PREREQ([2.2.2])
LT_INIT([disable-static dlopen])
AC_SUBST([LIBTOOL_DEPS])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

AC_LANG_PUSH([C++])

# Checks for libraries.
AC_ARG_WITH([blas-lapack],
            [AS_HELP_STRING([--with-blas-lapack=FLAGS],
            [Select the options to use to link the blas/lapack libraries.])],
            [BLAS_LIB=$withval],
            dnl look in /usr/lib/atlas
            [
                AC_MSG_NOTICE([searching for atlas in /usr/lib/atlas-base]); 
                save_LDFLAGS=$LDFLAGS; LDFLAGS="-L/usr/lib/atlas-base $LDFLAGS";
                AC_CHECK_LIB([f77blas],[daxpy_],[BLAS_LIB="-L/usr/lib/atlas-base -lf77blas -llapack_atlas -latlas"],
                dnl look in /usr/lib64/atlas
                [                                                          
                        AC_MSG_NOTICE([searching for atlas in /usr/lib64/atlas]);
                        LDFLAGS=$save_LDFLAGS; LDFLAGS="-L/usr/lib64/atlas $LDFLAGS";
                        AC_CHECK_LIB([f77blas],[daxpy_],[BLAS_LIB="-L/usr/lib64/atlas -lf77blas -llapack_atlas -latlas"],
                        dnl look in /opt/local/lib
                        [
                                AC_MSG_NOTICE([searching for atlas in /opt/local/lib]);
                                LDFLAGS=$save_LDFLAGS; LDFLAGS="-L/opt/local/lib $LDFLAGS";
                                AC_CHECK_LIB([f77blas],[daxpy_],[BLAS_LIB="-L/opt/local/lib -lf77blas -llapack -latlas"],cd cd ~
                                             [AC_MSG_ERROR([I can't find blas and lapack in standard places, please tell me where to search.])],
                                             [-latlas])
                        ], [-latlas])
                 ], [-latlas])
                 LDFLAGS=$save_LDFLAGS
            ])





AC_ARG_WITH([lislibrary],
						[AS_HELP_STRING([--without-lislibrary],
						[Disable lis library])],
						[],
						[with_lislibrary=yes])

AS_IF([test "x$with_lislibrary" != xno],
			[
				AC_ARG_WITH([lis-home],
						[AS_HELP_STRING([--with-lis-home=DIR],
						[Select the directory where LIS is installed.])],
						[LIS_HOME=$withval],
						[LIS_HOME=/opt/lis/1.4.60])

						LIS_CPPFLAGS="-DUSE_MPI -I$LIS_HOME/include"
						LIS_LDFLAGS="-L$LIS_HOME/lib"
						LIS_LIB="-llis"

						save_CPPFLAGS=$CPPFLAGS
						CPPFLAGS="$LIS_CPPFLAGS $CPPFLAGS"

						AC_CHECK_HEADER(lis.h, ,AC_MSG_ERROR([lis.h not found. I cowardly refuse to continue.]))

						save_LIBS=$LIBS
						save_LDFLAGS=$LDFLAGS
						LDFLAGS="$LIS_LDFLAGS $LDFLAGS"
						LIBS="$LIS_LIB $LIBS"
				AC_CHECK_LIB([lis],[lis_initialize],
					[AC_SUBST(LIS_CPPFLAGS, "$LIS_CPPFLAGS")
					 AC_SUBST(LIS_LDFLAGS, "$LIS_LDFLAGS")
					 AC_SUBST(LIS_LIB, "$LIS_LIB")
					 AC_DEFINE([HAVE_LIS],[1],[Define if you have lis library])
					 AC_MSG_NOTICE([It seems I am able to link LIS properly. Yuppee!!.])
					],
					[AC_MSG_ERROR([It seems I am not able to link LIS properly, --without-lislibrary to disable.])])])

AM_CONDITIONAL([HAVE_LIS],test x$with_lislibrary = xyes)

#AC_CHECK_LIB([lis], [lis_initialize], [AC_MSG_NOTICE([It seems I am able to link LIS properly. Yuppee!!.])],[AC_MSG_ERROR([It seems I am not able to link LIS properly. I cowardly refuse to continue.])])

LIBS=$save_LIBS
CPPFLAGS=$save_CPPFLAGS
LDFLAGS=$save_LDFLAGS

AC_ARG_WITH([mumps-home],
            [AS_HELP_STRING([--with-mumps-home=DIR],
            [Select the directory where MUMPS is installed.])],
            [MUMPS_HOME=$withval],
            [MUMPS_HOME=/usr])

MUMPS_CPPFLAGS="-I$MUMPS_HOME/include "
MUMPS_LDFLAGS="-L$MUMPS_HOME/lib -L$MUMPS_HOME/PORD/lib"
MUMPS_LIBS="-ldmumps -lmumps_common -lpord"

AC_ARG_WITH([mumps-extra-libs],
            [AS_HELP_STRING([--with-mumps-extra-libs=DIR],
            [Select extra libraries to link when linking to MUMPS.])],
            [MUMPS_EXTRA_LIBS=$withval],
            [MUMPS_EXTRA_LIBS=" -lptesmumps -lesmumps \
             -lscotch -lptscotch -lscotcherr -lptscotcherr -lscalapack-openmpi $BLAS_LIB \
              -lmpi_f90 -lmpi_f77 -lmpi -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.7 -lgfortran -lgfortranbegin"])


save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$MUMPS_CPPFLAGS $CPPFLAGS"
AC_CHECK_HEADER(dmumps_c.h, ,AC_MSG_ERROR([dmumps_c.h not found. I cowardly refuse to continue.]))

save_LIBS=$LIBS
save_LDFLAGS=$LDFLAGS
LDFLAGS="$MUMPS_LDFLAGS $LDFLAGS"
LIBS="$MUMPS_EXTRA_LIBS $LIBS"
 
AC_CHECK_LIB([dmumps], [dmumps_c], [AC_MSG_NOTICE([It seems I am able to link MUMPS properly. Yuppee!!.])],[AC_MSG_ERROR([It seems I am not able to link MUMPS properly. I cowardly refuse to continue.])],[-lmumps_common -lpord $MUMPS_EXTRA_LIBS])

LIBS=$save_LIBS
CPPFLAGS=$save_CPPFLAGS
LDFLAGS=$save_LDFLAGS

AC_SUBST(MUMPS_CPPFLAGS, "$MUMPS_CPPFLAGS")
AC_SUBST(MUMPS_LDFLAGS, "$MUMPS_LDFLAGS")
AC_SUBST(MUMPS_LIBS, "$MUMPS_LIBS")
AC_SUBST(MUMPS_EXTRA_LIBS, "$MUMPS_EXTRA_LIBS")
AC_SUBST(BLAS_LIB, "$BLAS_LIB")
#AC_SUBST(LIS_CPPFLAGS, "$LIS_CPPFLAGS")
#AC_SUBST(LIS_LDFLAGS, "$LIS_LDFLAGS")
#AC_SUBST(LIS_LIB, "$LIS_LIB")

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
