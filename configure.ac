# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.50])

AC_INIT([opmpressure], [0.1], [Bard.Skaflestad@sintef.no])

AM_INIT_AUTOMAKE([foreign -Wall])

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
AC_CONFIG_SRCDIR([src/mimetic.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h stddef.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memmove memset])

AX_LAPACK

AC_CONFIG_FILES([
  Makefile
  src/Makefile
])

AC_OUTPUT
