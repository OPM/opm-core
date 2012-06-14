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

OPM_CORE_CHECKS

OPM_DYNLINK_BOOST_TEST

AC_CONFIG_FILES([
 Makefile
 tests/Makefile
 examples/Makefile
 tutorials/Makefile
])

AC_OUTPUT
