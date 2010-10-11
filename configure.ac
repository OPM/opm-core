# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.50])

AC_INIT([libmimetic], [0.1], [Bard.Skaflestad@sintef.no])

AM_INIT_AUTOMAKE([foreign])

LT_INIT

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([mimetic.h])
AC_CONFIG_HEADERS([config.h])

# Neat trick.  Do we want it on by default?
AM_SILENT_RULES

# Checks for programs.
AC_PROG_CC
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
])

AC_OUTPUT
