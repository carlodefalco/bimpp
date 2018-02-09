# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_CONFIG_MACRO_DIR([m4])
AC_INIT([libbim], [0.0.1], [])
AC_CONFIG_SRCDIR([src/mesh.cpp])
AC_CONFIG_HEADERS([src/bim_config.h])
AH_TOP([#include "bim_config_pre.h"])

AM_INIT_AUTOMAKE([subdir-objects])

# Init libtool.
LT_PREREQ([2.2.2])
LT_INIT([disable-static dlopen])
AC_SUBST([LIBTOOL_DEPS])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

AC_LANG_PUSH([C++])

# Checks for libraries.
AC_F77_WRAPPERS
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

AC_ARG_WITH([p4estlibrary],
            [AS_HELP_STRING([--without-p4estlibrary],
            [Disable p4est library])],
            [],
            [with_p4estlibrary=yes])

AS_IF([test "x$with_p4estlibrary" != xno],
      [
       AC_ARG_WITH([p4est-home],
        [AS_HELP_STRING([--with-p4est-home=DIR],
        [Select the directory where p4est is installed.])],
        [P4EST_HOME=$withval],
        [P4EST_HOME=/usr/local])
      
        P4EST_CPPFLAGS="-DUSE_MPI -I$P4EST_HOME/include"
        P4EST_LDFLAGS="-L$P4EST_HOME/lib"
        P4EST_LIB="-lp4est -lsc"
        save_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$P4EST_CPPFLAGS $CPPFLAGS"

        AC_CHECK_HEADER(p4est.h, ,AC_MSG_ERROR([p4est.h not found. I cowardly refuse to continue.]))

        save_LIBS=$LIBS
        save_LDFLAGS=$LDFLAGS
        LDFLAGS="$P4EST_LDFLAGS $LDFLAGS"
        LIBS="$P4EST_LIB $LIBS"
        AC_CHECK_LIB([p4est],[p4est_refine],
        [AC_SUBST(P4EST_CPPFLAGS, "$P4EST_CPPFLAGS")
        AC_SUBST(P4EST_LDFLAGS, "$P4EST_LDFLAGS")
        AC_SUBST(P4EST_LIB,"$P4EST_LIB")
        AC_DEFINE([HAVE_P4EST],[1],[Define if you have p4est library])
        AC_MSG_NOTICE([It seems I am able to link p4est properly. Yuppee!!.])
       ],
       [AC_MSG_ERROR([It seems I am not able to link p4est properly. --without-p4estlibrary to disable])])])

AM_CONDITIONAL([HAVE_P4EST], test x$with_p4estlibrary = xyes)

LIBS=$save_LIBS
CPPFLAGS=$save_CPPFLAGS
LDFLAGS=$save_LDFLAGS


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
              -lmpi_f90 -lmpi_f77 -lmpi -lm"])


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


AC_ARG_WITH([octave-home],
            [AS_HELP_STRING([--with-octave-home=DIR],
            [Select the directory where octave mkoctfile is contained.])],
            [OCTAVE_HOME=$withval],
            [OCTAVE_HOME=/opt/octave/4.2.0/bin])

OCTAVE_BINDIR="$OCTAVE_HOME"

AC_PATH_PROG([OCTAVE_MKOCTFILE], [mkoctfile], [], ["${OCTAVE_BINDIR}"])

if [ test "${OCTAVE_MKOCTFILE}" != ""]; then


# use mkoctfile to get info on octave flags
   AC_SUBST(OCTAVE_COMPILE_FLAGS, "$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p INCFLAGS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p CPPFLAGS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p CXXPICFLAG) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p ALL_CXXFLAGS)")
   AC_SUBST(OCTAVE_LINK_FLAGS, "$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p ALL_CXXFLAGS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p RDYNAMIC_FLAG) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p ALL_LDFLAGS)")
   AC_SUBST(OCTAVE_LINK_ADD, "$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p LFLAGS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p OCTAVE_LIBS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p LAPACK_LIBS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p BLAS_LIBS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p FFTW_LIBS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p LIBREADLINE) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p LIBS) \
$(env -i PATH=$PATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH ${OCTAVE_MKOCTFILE} -p FLIBS) \
")
else
   AC_MSG_WARN([You do not seem to have a complete install of Octave on your computer.])
fi





# Checks for library functions.

AC_CONFIG_FILES([src/bim_config_pre.h
                 Makefile
                 src/Makefile
                 data/Makefile
                 test/Makefile
                 examples/tmesh_2d/Makefile
                 examples/adaptive_adr_2d/Makefile
                 examples/bgs/Makefile])
AC_OUTPUT

