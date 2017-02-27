# Process this script with autoconf to generate the configure script

dnl Standard initialization stuff
AC_INIT([overturns2d],[3.0],[vinodkl@stanford.edu],[overturns2d])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([subdir-objects])

dnl Fix installation directory to be source code directory
AC_PREFIX_DEFAULT('$(top_srcdir)')

dnl Setup standard programs
AC_PROG_AWK
AC_PROG_EGREP
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PROG_CC
AM_PROG_CC_C_O

dnl Change the order of the compiler selection
m4_define([MY_FC],["$FC" ifort gfortran g95])

AC_PROG_FC([MY_FC])

dnl Define F77 compiler and setting it to F90 compiler
AC_DEFINE(F77)
AC_SUBST(F77,["$FC"])

dnl declarations using appropriate compiler flags.
AC_ARG_ENABLE(real-double,
              AC_HELP_STRING([--enable-real-double],
                             [Use double precision as default for real (default: yes)]),,enable_real_double=yes)

AC_ARG_ENABLE(adjoints,
              AC_HELP_STRING([--enable-adjoints],
                             [Enable compilation of adjoints (default: yes)]),,enable_adjoints=yes)
AM_CONDITIONAL(COMPILE_ADJOINTS,[test x$enable_adjoints = xyes])

AC_ARG_ENABLE(utils,
              AC_HELP_STRING([--enable-utils],
                             [Enable compilation of utilities (default: yes)]),,enable_utils=yes)
AM_CONDITIONAL(COMPILE_UTILS,[test x$enable_utils = xyes])

AC_ARG_ENABLE(openmp,
              AC_HELP_STRING([--enable-openmp],
                             [Enable compilation with openmp (default: no)]),,enable_openmp=no)

dnl Set some useful compiler flags. However, if user has specified them 
dnl via command line, then do nothing
if test "$ac_cv_env_FCFLAGS_set" != set ; then
        FCFLAGS="-fPIC -O3"
fi

if test "$ac_cv_env_CFLAGS_set" != set ; then
        CFLAGS="-fPIC -O3"
fi

dnl Set up compiler flags for FC based on compiler chosen
case "$FC" in
        *ifort | *ifc)
                FCFLAGS="${FCFLAGS} -extend_source"
                if test x$enable_real_double = xyes ; then
                        FCFLAGS="${FCFLAGS} -r8"
                fi
                if test x$enable_openmp = xyes ; then
                        FCFLAGS="${FCFLAGS} -openmp"
                fi
                ;;
        *gfortran)
                FCFLAGS="${FCFLAGS} -ffixed-line-length-132"
                if test x$enable_real_double = xyes ; then
                        FCFLAGS="${FCFLAGS} -fdefault-real-8"
                fi
                if test x$enable_openmp = xyes ; then
                        FCFLAGS="${FCFLAGS} -fopenmp"
                fi
                ;;
        *g95)
                FCFLAGS="${FCFLAGS} -ffixed-line-length-132 -fno-second-underscore"
                if test x$enable_real_double = xyes ; then
                        FCFLAGS="${FCFLAGS} -r8"
                fi
                if test x$enable_openmp = xyes ; then
                        FCFLAGS="${FCFLAGS} -fopenmp"
                fi
                ;;
        *)
                AC_MSG_WARN([No preset configuration flags available for $FC. Set FCFLAGS such that default real is set to double precision and the 72 character limit is removed for .f files])
                ;;
esac

dnl Define F77 compiler flags and setting it to F90 flags
AC_DEFINE(FFLAGS)
AC_SUBST(FFLAGS,["$FCFLAGS"])

dnl Set up compiler flags for CC based on compiler chosen
case "$CC" in
        *icc)
                if test x$enable_openmp = xyes ; then
                        CFLAGS="$FCFLAGS} -openmp"
                fi
                ;;
        *gcc)
                if test x$enable_openmp = xyes ; then
                        CFLAGS="${CFLAGS} -fopenmp"
                fi
                ;;
        *g++)
                if test x$enable_openmp = xyes ; then
                        CFLAGS="${CFLAGS} -fopenmp"
                fi
                ;;
        *)
                AC_MSG_WARN([No preset configuration flags available for $CC.])
                ;;
esac

AC_MSG_NOTICE([Fortran compiler set to $FC])
AC_MSG_NOTICE([FCFLAGS set to $FCFLAGS])
AC_MSG_NOTICE([C compiler set to $CC])
AC_MSG_NOTICE([CFLAGS set to $FCFLAGS])

AC_LANG([Fortran])
AC_LANG_C

AC_CONFIG_FILES([Makefile
                 ihc/Makefile
                 ihc/src/Makefile
                 src/Makefile
                 src/adjoints/Makefile
                 utils/interpolation/Makefile
                 ])
AC_OUTPUT
