# whether compiler accepts -std=c++11 or -std=c++0x
# can be disabled by --disable-gxx0xcheck

AC_DEFUN([GXX0X],[
  ac_save_CXX="$CXX"
  
  # try flag -std=c++11
  AC_CACHE_CHECK([whether $CXX accepts -std=c++11], dune_cv_gplusplus_accepts_cplusplus11, [
    AC_REQUIRE([AC_PROG_CXX])
    AC_ARG_ENABLE(gxx0xcheck,
      AC_HELP_STRING([--disable-gxx0xcheck],
        [try flag -std=c++11 to enable C++11 features [[default=yes]]]),
        [gxx0xcheck=$enableval],
        [gxx0xcheck=yes])
    if test "x$GXX" = xyes && test "x$gxx0xcheck" = xyes; then
      AC_LANG_PUSH([C++])
      CXX="$CXX -std=c++11"
      AC_TRY_COMPILE([
        #include <iostream>
        #include <array>
        ], [],
        dune_cv_gplusplus_accepts_cplusplus11=yes,
        dune_cv_gplusplus_accepts_cplusplus11=no)
      AC_LANG_POP([C++])
    fi
  ])
  if test "x$dune_cv_gplusplus_accepts_cplusplus11" == "xyes" ; then
    CXX="$ac_save_CXX -std=c++11"
    CXXCPP="$CXXCPP -std=c++11"
  else
    CXX="$ac_save_CXX"
    
    # try flag -std=c++0x instead
    AC_CACHE_CHECK([whether $CXX accepts -std=c++0x], dune_cv_gplusplus_accepts_cplusplus0x, [
      AC_REQUIRE([AC_PROG_CXX])
      AC_ARG_ENABLE(gxx0xcheck,
        AC_HELP_STRING([--disable-gxx0xcheck],
          [try flag -std=c++0x to enable C++11 features [[default=yes]]]),
          [gxx0xcheck=$enableval],
          [gxx0xcheck=yes])
      if test "x$GXX" = xyes && test "x$gxx0xcheck" = xyes; then
        AC_LANG_PUSH([C++])
        CXX="$CXX -std=c++0x"
        AC_TRY_COMPILE([
          #include <iostream>
          #include <array>
          ], [],
          dune_cv_gplusplus_accepts_cplusplus0x=yes,
          dune_cv_gplusplus_accepts_cplusplus0x=no)
        AC_LANG_POP([C++])
      fi
    ])
    if test "x$dune_cv_gplusplus_accepts_cplusplus0x" == "xyes" ; then
      CXX="$ac_save_CXX -std=c++0x"
      CXXCPP="$CXXCPP -std=c++0x"
    else
      CXX="$ac_save_CXX"
    fi
  fi
])
