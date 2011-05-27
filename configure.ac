#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([cpgpreprocess], [1.0], [Jostein.R.Natvig@sintef.no])
AC_CONFIG_SRCDIR([preprocess.c])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([foreign])

# Use "silent" rules by default if available in the Automake being used to
# process this setup.  See, e.g.,
#
#   http://www.flameeyes.eu/autotools-mythbuster/automake/silent.html
#
# for details, particularly concerning the backwards compatibility.
#
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])


LT_INIT

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
