# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.50])

AC_INIT([opmtransport], [0.1], [Jostein.R.Natvig@sintef.no])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Use "silent" rules by default if available in the Automake being used to
# process this setup.  See, e.g.,
#
#   http://www.flameeyes.eu/autotools-mythbuster/automake/silent.html
#
# for details, particularly concerning the backwards compatibility.
#
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

LT_INIT

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/spu_explicit.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX

AC_PROG_F77
AC_F77_WRAPPERS

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AX_LAPACK

# Checks for header files.
AC_CHECK_HEADERS([assert.h limits.h stdlib.h])

AC_CHECK_HEADERS([sparse_sys.h],
                 [sparse_sys_header=yes],
                 [sparse_sys_header=no])

AC_CHECK_HEADERS([suitesparse/umfpack.h],
                 [umfpack_header=yes],
                 [umfpack_header=no])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memmove memset pow sqrt])

AC_SEARCH_LIBS([umfpack_dl_solve], [umfpack],
               [umfpack_lib=yes], [umfpack_lib=no])

AC_SEARCH_LIBS([csrmatrix_elm_index], [opmpressure],
               [sparse_sys_lib=yes],
               [sparse_sys_lib=no], [$LAPACK_LIBS $BLAS_LIBS])

AM_CONDITIONAL([UMFPACK],
  [test "x$umfpack_header" != "xno" -a "x$umfpack_lib" != "xno"])
AM_CONDITIONAL([OPMPRESSURE_SPARSE_SYS],
  [test "x$sparse_sys_header" != "xno" -a "x$sparse_sys_lib" != "xno"])

AC_CONFIG_FILES([
  Makefile
  src/Makefile
])

AC_OUTPUT
