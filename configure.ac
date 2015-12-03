#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([GBMUnTrig], [1.2], [chapliv@uah.edu], [GBMUnTrig])

AM_INIT_AUTOMAKE([GBMUnTrig], [1.2])

AC_PREFIX_DEFAULT([/usr/local/GBMUnTrig])

AC_CONFIG_SRCDIR([src/gbmCont2UnTrig.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for library functions.
AC_HEADER_STDC
AC_HEADER_STDBOOL
AC_CHECK_FUNCS([strptime])
  
AC_ARG_WITH([heasoft],
  [AS_HELP_STRING([--with-heasoft],
    [location of HEASOFT installation. Overwrites HEADAS. Assumes this directory contains lib/ and include/])],
  [LHEAS_CFLAGS="-I$withval/include"
   LHEAS_LDFLAGS="-L$withval/lib"
   noLHEASflags=false],
  [noLHEASflags=true])
 
AC_ARG_WITH([cfitsio],
  [AS_HELP_STRING([--with-cfitsio],
    [location of CFITSIO installation. This has the same effect as --with-heasoft.])],
  [LHEAS_CFLAGS="-I$withval/include"
   LHEAS_LDFLAGS="-L$withval/lib"
   noLHEASflags=false],
  [])
  
AC_ARG_WITH([cfitsio-include],
  [AS_HELP_STRING([--with-cfitsio-include],
    [location of fitsio.h if non-standard installation])],
  [LHEAS_CFLAGS="-I$withval/include"
   noLHEASflags=false],
  [])

AC_ARG_WITH([cfitsio-lib],
  [AS_HELP_STRING([--with-cfitsio-lib],
    [location of cfitsio library if non-standard installation])],
  [LHEAS_LDFLAGS="-L$withval/lib"
   noLHEASflags=false],
  [])

AC_CHECK_PROG(WHICH_RM, rm, rm)
AC_CHECK_PROG(WHICH_ENV, env, env)


AM_CONDITIONAL(PUTVERSN, test x$putvernum = xtrue)
AM_CONDITIONAL(USE_HEADAS, test x$noLHEASflags = xtrue)

AC_SUBST([LHEAS_CFLAGS])
AC_SUBST([LHEAS_LDFLAGS])

AC_ARG_VAR(HEADAS, Location of heasoft. Overwritten by '--with-heasoft'.)

AC_SUBST([HEADAS])


AC_PROG_MAKE_SET
AC_CONFIG_FILES([Makefile])
AC_OUTPUT