#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([OPM Core Library], [0.1], [atgeirr@sintef.no],
        [opmcore], [https://public.ict.sintef.no/opm/hg/opmcore])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Needed for automake since version 1.12 because extra-portability
# warnings were then added to -Wall. Ifdef makes it backwards compatible.
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([opm/core/grid.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

dnl Initialize libtool; the funny indentation here is to
dnl satisfy libtoolize' check for the presence of this macro
m4_ifdef([LT_INIT],
         [
          LT_INIT[]dnl
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

dnl Substitute Autoconf's abs_*dir variables into the Makefiles for the
dnl benefit of external code that uses these variables to derive
dnl locations (e.g., Dune's DUNE_CHECK_MODULES macro).  Automakes prior
dnl to version 1.10 do not automatically substitute these variables into
dnl output files.
AC_SUBST([abs_srcdir])
AC_SUBST([abs_builddir])
AC_SUBST([abs_top_srcdir])
AC_SUBST([abs_top_builddir])

AC_CONFIG_FILES([
 Makefile
 tests/Makefile
 examples/Makefile
 tutorials/Makefile
 opm-core.pc
 lib/pkgconfig/opm-core.pc
])

AC_OUTPUT
