dnl #################################################################
dnl File: ./configure.ac
dnl Package: landsepi
dnl Description:
dnl   + autoconf file
dnl
dnl Process this file with autoconf to produce a configure script.
dnl #################################################################

dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------

define([THE_PACKAGE_NAME],    [landsepi])
define([THE_PACKAGE_VERSION], [1.5.3])
define([THE_PACKAGE_MAINTAINER],  [jean-francois.rey@inrae.fr])

dnl --------------------------------
dnl Initialization macros.
dnl --------------------------------
AC_INIT([THE_PACKAGE_NAME],[THE_PACKAGE_VERSION],[THE_PACKAGE_MAINTAINER])

dnl -----------------------------------------------
dnl Checks for programs.
dnl -----------------------------------------------

dnl check R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR(["could not determine R_HOME, is R installed ?"])
fi
AC_MSG_NOTICE([${R_HOME}])
#CXX1X=`"${R_HOME}/bin/R" CMD config CXX1X`
#CXX1XSTD=`"${R_HOME}/bin/R" CMD config CXX1XSTD`
#CXX="${CXX1X} ${CXX1XSTD}"
CXX=`"${R_HOME}/bin/R" CMD config CXX`
#CPP=`"${R_HOME}/bin/R" CMD config CPP`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
#CXX1XFLAGS=`"${R_HOME}/bin/R" CMD config CXX1XFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

AC_LANG(C++)

## debug mode for compiler code 
AC_MSG_CHECKING([whether to build debug version])
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[build debug version (default=no)]))

if test "$enable_debug" != "yes" ; then
  AC_MSG_RESULT(no)
  CXXFLAGS="${CXXFLAGS}"
else
  AC_MSG_RESULT(yes)
  CXXFLAGS="${CXXFLAGS} -pedantic -Wall -Wconversion -Wextra -g -ggdb -O0 -fno-inline"
  CPPFLAGS=" -UNDEBUG -DDEBUG=1 ${CPPFLAGS}"
fi

## profiling mode for compiler code 
AC_MSG_CHECKING([whether to build profiling version])
AC_ARG_ENABLE(prof, AS_HELP_STRING([--enable-prof],[build profiling version (default=no)]))

if test "$enable_prof" != "yes" ; then
  AC_MSG_RESULT(no)
else
  AC_MSG_RESULT(yes)
  #AC_CHECK_LIB([profiler])
  AC_PATH_PROG([Gprof], [google-pprof])

  if test "${Gprof}" != ""; then
    GPROF_LIBS="-lprofiler"
  else
    AC_MSG_ERROR([google-perftools/gperftools not found, is google perftools installed?])
  fi

fi



# check Rcpp
AC_MSG_NOTICE([${R_HOME}])
AC_MSG_CHECKING([for Rcpp package])
PKG_RCPP_LIBS=`${R_HOME}/bin/Rscript -e "Rcpp:::LdFlags()"`
PKG_RCPP_CXXFLAGS=`${R_HOME}/bin/Rscript -e "Rcpp:::CxxFlags()"`
AC_MSG_RESULT(yes)
AC_SUBST(PKG_RCPP_LIBS)
AC_SUBST(PKG_RCPP_CXXFLAGS)

# GSL librarie 
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

LDFLAGS="${LDFLAGS} ${LIBS} ${GPROF_LIBS}"
CXXFLAGS="${CXXFLAGS}"


dnl -----------------------------------------------
dnl Generates Makefile's, configuration files and scripts
dnl -----------------------------------------------
# R mode
AC_SUBST(CXX)
#AC_SUBST(CXX11)
AC_SUBST(CXXFLAGS)
#AC_SUBST(CXX1XFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
