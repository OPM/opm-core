AC_DEFUN([OPM_LAPACK],
[AC_REQUIRE([AC_F77_WRAPPERS])dnl
 AC_REQUIRE([AX_LAPACK])dnl
 if test x"$ax_lapack_ok" != xyes; then
 	AC_MSG_ERROR([BLAS/LAPACK required, but not found.])
 fi
])[]dnl
