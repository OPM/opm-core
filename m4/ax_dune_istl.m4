AC_DEFUN([OPM_DUNE_ISTL_PROGRAM_TEXT],
[AC_LANG_PROGRAM(
  [[#include <dune/common/deprecated.hh>
    #include <dune/istl/bcrsmatrix.hh>
    #include <dune/common/fmatrix.hh>
  ]],dnl
  [[typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > Matrix;
    Matrix matrix( 3, 3, Matrix::random );
    for (int i = 0; i < 3; ++i) matrix.setrowsize(i, 2);
    matrix.endrowsizes();
  ]])[]dnl
])


AC_DEFUN([AX_DUNE_ISTL],
[AC_REQUIRE([AX_DUNE_COMMON])

 AC_CACHE_CHECK(dnl
[for installed dune-istl headers],dnl
[ax_cv_dune_istl_available],dnl
  [AC_LANG_PUSH([C++])[]dnl

   AC_LINK_IFELSE([OPM_DUNE_ISTL_PROGRAM_TEXT],dnl
                  [ax_cv_dune_istl_available=yes],dnl
                  [ax_cv_dune_istl_available=no])

   AC_LANG_POP([C++])[]dnl
  ])[]dnl

  AS_IF([test "x$ax_cv_dune_istl_available" = "xyes" -a \
              "x$ax_cv_dune_common_available" = "xyes"],dnl
   [AC_DEFINE([HAVE_DUNE_ISTL], [1],dnl
              [Define to 1 if `dune-istl' is available])
   ])[]dnl

  AM_CONDITIONAL([DUNE_ISTL],
    [test "x$ax_cv_dune_istl_available" = "xyes" -a \
          "x$ax_cv_dune_common_available" = "xyes"])
])
