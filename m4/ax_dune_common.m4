AC_DEFUN([OPM_DUNE_COMMON_PROGRAM_TEXT],
[AC_LANG_PROGRAM(
  [[#include <dune/common/fvector.hh>
    #include <dune/common/fmatrix.hh>
  ]],dnl
  [[Dune::FieldVector<double,1>   v;
    Dune::FieldMatrix<double,1,1> m;
    m[0][0] = 1.0;
    v[0]    = 1.0;
    Dune::FieldVector<double,1> w = m*v;
  ]])[]dnl
])


AC_DEFUN([AX_DUNE_COMMON],
[AC_CACHE_CHECK(dnl
[for installed dune-common headers],dnl
[ax_cv_dune_common_available],dnl
  [AC_LANG_PUSH([C++])[]dnl

   AC_LINK_IFELSE([OPM_DUNE_COMMON_PROGRAM_TEXT],dnl
                  [ax_cv_dune_common_available=yes],dnl
                  [ax_cv_dune_common_available=no])

   AC_LANG_POP([C++])[]dnl
  ])[]dnl

  AS_IF([test "x$ax_cv_dune_common_available" = "xyes"],dnl
   [AC_DEFINE([HAVE_DUNE_COMMON], [1],dnl
              [Define to 1 if `dune-common' is available])[]dnl
    LIBS="-ldunecommon $LIBS"
   ])[]dnl

  AM_CONDITIONAL([DUNE_COMMON],
   [test "x$ax_cv_dune_common_available" = "xyes"])
])
