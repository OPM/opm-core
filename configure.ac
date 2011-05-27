# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.50])

AC_INIT([opmtransport], [0.1], [Jostein.R.Natvig@sintef.no])

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
AC_CONFIG_SRCDIR([src/spu_explicit.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([assert.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memmove memset])

AC_CONFIG_FILES([
  Makefile
  src/Makefile
])

AC_OUTPUT
