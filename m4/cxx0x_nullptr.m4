AC_DEFUN([NULLPTR_CHECK],[
  AC_CACHE_CHECK([whether nullptr is supported], dune_cv_nullptr_support, [
    AC_REQUIRE([AC_PROG_CXX])
    AC_REQUIRE([GXX0X])
    AC_LANG_PUSH([C++])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
      char* ch = nullptr;
      if(ch!=nullptr) { ; }
      ])],
      [dune_cv_nullptr_support=yes],
      [dune_cv_nullptr_support=no])
    AC_LANG_POP
  ])
  if test "x$dune_cv_nullptr_support" = xyes; then
    AC_DEFINE(HAVE_NULLPTR, 1, [Define to 1 if nullptr is supported])
  fi
])
