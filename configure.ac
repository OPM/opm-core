#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.62])
AC_INIT([OPM Core Library], [0.1], [atgeirr@sintef.no],dnl
        [opmcore], [https://public.ict.sintef.no/opm/hg/opmcore])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

LT_INIT

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([opmcore/grid.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CC_C_O
AC_PROG_LIBTOOL

# AX_LAPACK requires a working F77 compiler or, rather, its runtime
# support libraries.
AC_PROG_F77

# F77 name mangling
AC_F77_WRAPPERS

# Checks for libraries.
AX_LAPACK

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
AC_CHECK_FUNCS([floor memmove pow sqrt strchr strtol])

AC_SEARCH_LIBS([umfpack_dl_solve], [umfpack],
               [umfpack_lib=yes], [umfpack_lib=no])

AM_CONDITIONAL([UMFPACK],
  [test "x$umfpack_header" != "xno" -a "x$umfpack_lib" != "xno"])


AC_CONFIG_FILES([
 Makefile
 tests/Makefile
 examples/Makefile
])

AC_OUTPUT
