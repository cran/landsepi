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
define([THE_PACKAGE_VERSION], [0.0.2])
define([THE_PACKAGE_MAINTAINER],  [jean-francois.rey@inra.fr])

dnl --------------------------------
dnl Initialization macros.
dnl --------------------------------
AC_INIT(THE_PACKAGE_NAME, THE_PACKAGE_VERSION, THE_PACKAGE_MAINTAINER)

dnl -----------------------------------------------
dnl Checks for programs.
dnl -----------------------------------------------

dnl check R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR(["could not determine R_HOME, is R installed ?"])
fi
AC_MSG_NOTICE([${R_HOME}])
CXX1X=`"${R_HOME}/bin/R" CMD config CXX1X`
CXX1XSTD=`"${R_HOME}/bin/R" CMD config CXX1XSTD`
CXX="${CXX1X} ${CXX1XSTD}"
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXX1XFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

AC_LANG(C++)

## debug mode for compiler code 
AC_MSG_CHECKING([whether to build debug version])
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug], [build debug version (default=no)]))

if test "$enable_debug" != "yes" && test "$DEBUG" != TRUE ; then
  AC_MSG_RESULT(no)
  CXXFLAGS="${CXXFLAGS} -w"
else
  AC_MSG_RESULT(yes)
  CXXFLAGS="${CXXFLAGS} -Wall -Wextra -g -ggdb -Wconversion -O0 -fno-inline -DDEBUG"
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

LDFLAGS="${LDFLAGS} ${LIBS}"

# GDAL librarie
GDAL_CONFIG="gdal-config"

AC_MSG_CHECKING(gdal-config usability)
if test `${GDAL_CONFIG} --version`;
then
  GDAL_CFLAGS=`${GDAL_CONFIG} --cflags`
  GDAL_LIBS=`${GDAL_CONFIG} --libs`
  GDAL_VERSION=`${GDAL_CONFIG} --version`
  GDAL_DEP_LIBS=`${GDAL_CONFIG} --dep-libs`
  GDAL_DATADIR=`${GDAL_CONFIG} --datadir`
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  echo "Error: gdal-config not found"
  echo "The gdal-config script distributed with GDAL could not be found."
  echo "If you have not installed the GDAL libraries, you can"
  echo "download the source from  http://www.gdal.org/"
  echo "If you have installed the GDAL libraries, then make sure that"
  echo "gdal-config is in your path. Try typing gdal-config at a"
  echo "shell prompt and see if it runs. If not, use:"
  echo " --configure-args='--with-gdal-config=/usr/local/bin/gdal-config'"
  echo "with appropriate values for your installation."
  echo ""
  exit 1

fi

AC_MSG_NOTICE([GDAL: ${GDAL_VERSION}])
AC_MSG_CHECKING([GDAL version >= 1.11.0])
GDAL_VER_MAX=`echo $GDAL_VERSION | cut -d'.' -f1`
GDAL_VER_MED=`echo $GDAL_VERSION | cut -d'.' -f1`
if test ${GDAL_VER_MAX} -lt 1 && test ${GDAL_VER_MED} -lt 11; then
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([upgrade GDAL to 1.11.0 or later])
else
  AC_MSG_RESULT(yes)
  CPPFLAGS="${CPPFLAGS} -DGDAL_VERSION=${GDAL_VER_MAX}"
  CXXFLAGS="${CXXFLAGS} ${GDAL_CFLAGS}"
fi

LDFLAGS="${LDFLAGS} ${GDAL_LIBS}"

CXXFLAGS="${CXXFLAGS}"

dnl -----------------------------------------------
dnl Generates Makefile's, configuration files and scripts
dnl -----------------------------------------------
# R mode
AC_SUBST(CXX)
AC_SUBST(CXXFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
