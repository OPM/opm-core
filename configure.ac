#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([OPM Core Library], [0.1], [atgeirr@sintef.no],dnl
        [opmcore], [https://public.ict.sintef.no/opm/hg/opmcore])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

m4_ifdef([LT_INIT], [LT_INIT])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([opm/core/grid.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
m4_ifdef([LT_INIT], [], [AC_PROG_LIBTOOL])
AC_PROG_RANLIB

# AX_LAPACK requires a working F77 compiler or, rather, its runtime
# support libraries.
AC_PROG_F77

# F77 name mangling
AC_F77_WRAPPERS

# Checks for libraries.

# Bring in numerics support (standard library component)
AC_SEARCH_LIBS([sqrt], [m])

AX_LAPACK
AX_BOOST_BASE([1.37])
AX_BOOST_SYSTEM
AX_BOOST_DATE_TIME
AX_BOOST_FILESYSTEM
AX_BOOST_UNIT_TEST_FRAMEWORK

AX_DUNE_ISTL

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


AC_CONFIG_FILES([
 Makefile
 tests/Makefile
 examples/Makefile
 tutorials/Makefile
])

AC_OUTPUT
