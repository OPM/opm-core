AC_DEFUN([STATIC_ASSERT_CHECK],[
  AC_CACHE_CHECK([whether static_assert is supported], dune_cv_static_assert_support, [
    AC_REQUIRE([AC_PROG_CXX])
    AC_REQUIRE([GXX0X])
    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[static_assert(true,"MSG")])],
      [dune_cv_static_assert_support=yes],
      [dune_cv_static_assert_support=no])
    AC_LANG_POP
  ])
  if test "x$dune_cv_static_assert_support" = xyes; then
    AC_DEFINE(HAVE_STATIC_ASSERT, 1, [Define to 1 if static_assert is supported])
  fi
])
