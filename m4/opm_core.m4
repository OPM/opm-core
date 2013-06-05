dnl -*- autoconf -*-

dnl locate opm-core library itself; this macro is called by every module
dnl that depends on opm-core.
AC_DEFUN([OPM_CORE_CHECK_MODULE],
[
 OPM_CHECK_PKG_MODULE([opm-core],[1.0],[OPM Core Library])
])

dnl find all prerequisites of opm-core; nothing to do here since this
dnl is done by the CMake module and then stored in the -config file.
AC_DEFUN([OPM_CORE_CHECKS],[])
