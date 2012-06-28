#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([OPM Core Library], [0.1], [atgeirr@sintef.no],dnl
        [opmcore], [https://public.ict.sintef.no/opm/hg/opmcore])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([opm/core/grid.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

m4_ifdef([LT_INIT],
         [LT_INIT[]dnl
          LT_LANG([C++])dnl
          LT_LANG([Fortran 77])dnl
          LT_LANG([Fortran])dnl
         ],dnl
         [AC_PROG_LIBTOOL[]dnl
          AC_PROG_CXX[]dnl
          AC_PROG_F77[]dnl
          AC_PROG_FC[]dnl
         ])[]dnl

# Checks for libraries.

# Bring in numerics support (standard library component)
AC_SEARCH_LIBS([sqrt], [m])

OPM_LAPACK

AX_BOOST_BASE([1.37])
AX_BOOST_SYSTEM
AX_BOOST_DATE_TIME
AX_BOOST_FILESYSTEM
AX_BOOST_UNIT_TEST_FRAMEWORK

AX_DUNE_ISTL
OPM_AGMG

OPM_DYNLINK_BOOST_TEST

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stddef.h stdlib.h string.h])

AC_CHECK_HEADERS([suitesparse/umfpack.h],
                 [umfpack_header=yes],
                 [umfpack_header=no])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_CHECK_FUNCS([floor memset memmove strchr strtol sqrt pow])
AC_FUNC_STRTOD

# Search for UMFPACK direct sparse solver.
AC_SEARCH_LIBS([amd_free],             [amd])
AC_SEARCH_LIBS([camd_free],            [camd])
AC_SEARCH_LIBS([colamd_set_defaults],  [colamd])
AC_SEARCH_LIBS([ccolamd_set_defaults], [ccolamd])
AC_SEARCH_LIBS([cholmod_l_start],      [cholmod])
AC_SEARCH_LIBS([umfpack_dl_solve],     [umfpack],dnl
               [umfpack_lib=yes],      [umfpack_lib=no])

AM_CONDITIONAL([UMFPACK],
  [test "x$umfpack_header" != "xno" -a "x$umfpack_lib" != "xno"])

m4_ifdef([AM_COND_IF],
[AM_COND_IF([UMFPACK], [],
 [AC_MSG_NOTICE([Found no working installation of UMFPACK.
                 UMFPACK support is disabled.])])
])

#-----------------------------------------------------------------
# Support for ERT configuration in the OPM build process. The
# configuration is unfortunately not really functional :-(.
#
# Define the option --with-ert=/path/to/ert/root which can be passed
# to the configure script. If we pass --with-ert=no (the default)
# support for use of ert libraries will not be added built in;
# otherwise we will build support for ERT. 
#

# The argument given to the --with-ert option will be interpreted as
# the root path of a ert installation, and the build process will
# expect to find the ert libraries in $with_ert/lib and the header
# files in $with_ert/include. There are no check for actually finding
# the required libraries and headers, and support for ert could in
# principle be configured as:
#
# configure --with-ert=yes CPPFLAGS=-I/path/to/ert/include  LDFLAGS=-L/path/to/ert/lib

AC_ARG_WITH(ert, [AS_HELP_STRING([--with-ert=<root>], [Use ERT libraries])], [], [with_ert=no])

AM_CONDITIONAL([HAVE_ERT], 
   [test "x$ert" != "xno"])

# The purpose of this macro was to test if the WITH_ERT conditional
# was true, and in that case the statements:
#
#   AC_SUBST([ERT_ROOT],["$with_ert"])
#   AC_DEFINE([HAVE_ERT],1,[Are the ERT libraries available])
# 
# should be invoked - but the if construction does not work, and at
# the bottom of the configuration the two statements are invoked
# unconditionally (i.e. in the end the build will be with ert
# anyway). I will award anyone who can fix this a real bear :-) -
# Joakim.
m4_ifdef([AM_COND_IF],
 [AM_COND_IF([HAVE_ERT], [AC_SUBST([ERT_ROOT],["$with_ert"])], [])
])

AC_SUBST([ERT_ROOT],["$with_ert"])
AC_DEFINE(HAVE_ERT,1,[Are the ERT libraries available for reading and writing ECLIPSE files.])

# In addition to what is present in this script ERT support is further
# configured through:
#
#  o The ERT_LDFLAGS and ERT_CPPFLAGS configured in the Makefile.am
#    template.
#
#  o The config.h file will define the symbol HAVE_ERT if ert is built into
#    the current installation.
#-----------------------------------------------------------------



AC_CONFIG_FILES([
 Makefile
 tests/Makefile
 examples/Makefile
 tutorials/Makefile
])

AC_OUTPUT
