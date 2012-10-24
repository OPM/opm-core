# whether compiler accepts -std=c++11 or -std=c++0x
# can be disabled by --disable-gxx0xcheck

AC_DEFUN([GXX0X],[
  save_CXX="$CXX"

  # put this check first, so we get disable C++11 if C++0x is
  AC_ARG_ENABLE(gxx0xcheck,
    AC_HELP_STRING([--disable-gxx0xcheck],
      [try flag -std=c++0x to enable C++0x features [[default=yes]]]),
      [gxx0xcheck=$enableval],
      [gxx0xcheck=yes])
  AC_ARG_ENABLE(gxx11check,
    AC_HELP_STRING([--disable-gxx11check],
      [try flag -std=c++11 to enable C++11 features [[default=yes]]]),
      [gxx11check=$enableval],
      [gxx11check=yes])

  # try flag -std=c++11
  AC_CACHE_CHECK([whether $CXX accepts -std=c++11], dune_cv_gplusplus_accepts_cplusplus11, [
    AC_REQUIRE([AC_PROG_CXX])
    if test "x$GXX" = xyes && test "x$gxx11check" = xyes && test "x$gxx0xcheck" = xyes ; then
      AC_LANG_PUSH([C++])
      CXX="$CXX -std=c++11"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
        #include <iostream>
        #include <array>
        ],)],
        [dune_cv_gplusplus_accepts_cplusplus11=yes],
        [dune_cv_gplusplus_accepts_cplusplus11=no])
      AC_LANG_POP([C++])
    else
      dune_cv_gplusplus_accepts_cplusplus11=disabled
    fi
  ])
  if test "x$dune_cv_gplusplus_accepts_cplusplus11" == "xyes" ; then
    CXX="$save_CXX -std=c++11"
    CXXCPP="$CXXCPP -std=c++11"
  else
    CXX="$save_CXX"
    
    # try flag -std=c++0x instead
    AC_CACHE_CHECK([whether $CXX accepts -std=c++0x], dune_cv_gplusplus_accepts_cplusplus0x, [
      AC_REQUIRE([AC_PROG_CXX])
      if test "x$GXX" = xyes && test "x$gxx0xcheck" = xyes; then
        AC_LANG_PUSH([C++])
        CXX="$CXX -std=c++0x"
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
          #include <iostream>
          #include <array>
          ],)],
          [dune_cv_gplusplus_accepts_cplusplus0x=yes],
          [dune_cv_gplusplus_accepts_cplusplus0x=no])
        AC_LANG_POP([C++])
      else
        dune_cv_gplusplus_accepts_cplusplus0x=disabled
      fi
    ])
    if test "x$dune_cv_gplusplus_accepts_cplusplus0x" == "xyes" ; then
      CXX="$save_CXX -std=c++0x"
      CXXCPP="$CXXCPP -std=c++0x"
    else
      CXX="$save_CXX"
    fi
  fi
])
