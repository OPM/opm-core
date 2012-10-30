## -*- autoconf -*-
# $Id: superlu.m4 5908 2010-02-23 16:46:05Z joe $
# searches for SuperLU headers and libs

# _slu_lib_path(SUPERLU_ROOT, HEADER)
#
# Try to find the subpath unter SUPERLU_ROOT containing HEADER.  Try
# SUPERLU_ROOT/"include/superlu", SUPERLU_ROOT/"include", and
# SUPERLU_ROOT/"SRC", in that order.  Set the subpath for the library to
# "lib".  If HEADER was found in SUPERLU_ROOT/"SRC", check whether
# SUPERLU_ROOT/"lib" is a directory, and set the subpath for the library to
# the empty string "" if it isn't.
#
# Shell variables:
#   my_include_path
#      The subpath HEADER was found in: "include/superlu", "include", or
#      "SRC".  Contents is only meaningful for my_slu_found=yes.
#   my_lib_path
#      The subpath for the library: "lib" or "". Contents is only meaningful
#      for my_slu_found=yes.
#   my_slu_found
#      Whether HEADER was found at all.  Either "yes" or "no".
AC_DEFUN([_slu_lib_path],
    [
        my_include_path=include/superlu
        my_lib_path=lib
        my_slu_found=yes
        if test ! -f "$1/$my_include_path/$2" ; then
            #Try to find headers under superlu
            my_include_path=include
            if test ! -f "$1/$my_include_path/$2" ; then
                my_include_path=SRC
                if test ! -f "$1/$my_include_path/$2"; then
                    my_slu_found=no
                else
                    if ! test -d "$1/$my_lib_path"; then
                        my_lib_path=""
                    fi
                fi
            fi
        fi
    ]
)

# _slu_search_versions(SUPERLU_ROOT)
#
# Search for either "slu_ddefs.h" or "dsp_defs.h" using _slu_lib_path().
#
# Shell variables:
#   my_slu_header
#      The name of the header that was found: first of "slu_ddefs.h" or
#      "dsp_defs.h".  Contents is only meaningful for my_slu_found=yes.
#   my_include_path
#      The subpath the header was found in: "include/superlu", "include", or
#      "SRC".  Contents is only meaningful for my_slu_found=yes.
#   my_lib_path
#      The subpath for the library: "lib" or "". Contents is only meaningful
#      for my_slu_found=yes.
#   my_slu_found
#      Whether any of the headers.  Either "yes" or "no".
AC_DEFUN([_slu_search_versions],
    [
        my_slu_header=slu_ddefs.h
        _slu_lib_path($1, $my_slu_header)
        if test "$my_slu_found" != "yes"; then
            my_slu_header="dsp_defs.h"
            _slu_lib_path($1, $my_slu_header)
        fi
    ]
)


# _slu_search_default()
#
# Search for SuperLU in the default locations "/usr" and "/usr/local".
#
# Shell variables:
#   with_superlu
#     Root of the SuperLU installation: first of "/usr" and "/usr/local".
#     Contents is only meaningful for my_slu_found=yes.
#   For other output variables see documentation of _slu_search_versions().
AC_DEFUN([_slu_search_default],
    [
        with_superlu=/usr
        _slu_search_versions($with_superlu)

        if test "$my_slu_found" = "no"; then
            with_superlu=/usr/local
            _slu_search_versions($with_superlu)
        fi
    ]
)


# OPM_PATH_SUPERLU()
#
# REQUIRES: AC_PROG_CC, AX_BLAS
#
# Shell variables:
#   with_superlu
#     "no", "yes (version 4.3 or newer)", "yes (version 4.2 or older, post 2005)" or "yes (pre 2005)"
#   direct_SUPERLU_CPPFLAGS
#   direct_SUPERLU_LIBS
#      CPPFLAGS and LIBS necessary to link against SuperLU.  This variable
#      contains no indirect references and is suitable for use inside
#      configure.  Guaranteed empty if SuperLU was not found.
#   SUPERLU_CPPFLAGS
#   SUPERLU_LIBS
#      CPPFLAGS and LIBS necessary to link against SuperLU.  This variable may
#      contain indirect references and is suitable for use inside makefiles.
#      Guaranteed empty if SuperLU was not found.
#   HAVE_SUPERLU
#      "0" or "1" depending on whether SuperLU was found.
#
# Substitutions:
#   SUPERLU_LIBS
#   SUPERLU_CPPFLAGS
#     Substitutes the values of the corresponding shell variables.
#   ALL_PKG_LIBS
#   ALL_PKG_CPPFLAGS
#     Adds references to SuperLU's substitutions.
#
# Defines:
#   HAVE_SUPERLU
#     ENABLE_SUPERLU or undefined.  Whether SuperLU was found.  The correct
#     way to check this is "#if HAVE_SUPERLU": This way SuperLU features will
#     be disabled unless $(SUPERLU_CPPFLAGS) was given when compiling.
#   SUPERLU_POST_2005_VERSION
#     1 or undefined.  A post-2005 version of SuperLU uses the header
#     "slu_ddefs.h" while earlier versions use "dsp_defs.h".
#   SUPERLU_MIN_VERSION_4_3
#     1 or undefined. SuperLU version 4.3 or newer uses the symbol
#     "SLU_DOUBLE" while earlier versions use "DOUBLE".
#   HAVE_MEM_USAGE_T_EXPANSIONS
#     1 or undefined.  Whether "mem_usage_t.expansions" was found in
#     "slu_ddefs.h" or "dsp_defs.h" as apropriate.
#
# Conditionals:
#   SUPERLU
AC_DEFUN([OPM_PATH_SUPERLU],[
    AC_REQUIRE([AC_PROG_CC])
    # we need this for FLIBS
    AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
    AC_REQUIRE([AX_BLAS])

    #
    # User hints ...
    #
    my_lib_path=""
    my_include_path=""
    AC_ARG_WITH([superlu],
        [AC_HELP_STRING([--with-superlu],[user defined path to SuperLU library])],
        [dnl
            if test x"$withval" != xno ; then
                # get absolute path
                with_superlu=`eval cd $withval > /dev/null 2>&1 && pwd`
                if test x"$withval" = xyes; then
                    # Search in default locations
                    _slu_search_default
                else
                    # Search for the headers in the specified location
                    _slu_search_versions(["$with_superlu"])
                fi
            fi
        ], [dnl
            # Search in default locations
            _slu_search_default
        ])

    AC_ARG_WITH([superlu-lib],
        [AC_HELP_STRING([--with-superlu-lib],
                [The name of the static SuperLU library to link to.  By default
                    the static library with the name superlu.a is tried, but
                    only if shared linking has failed first.  Giving this
                    options forces static linking for SuperLU.])],
        [
            if test x"$withval" = xno ; then
                with_superlu_lib=
            fi
        ])

    AC_ARG_WITH([superlu-blaslib],
        [AC_HELP_STRING([--with-superlu-blaslib],
                [The name of the static blas library to link to.  By default
                    we try to link to the systems blas library (see
                    --with-blas).  Giving this options forces static linking
                    for SuperLU.])],
        [
            if test "$withval" = no ; then
                with_superlu_blaslib=
            fi
        ])

    # store old values
    ac_save_LDFLAGS="$LDFLAGS"
    ac_save_CPPFLAGS="$CPPFLAGS"
    ac_save_LIBS="$LIBS"

    # do nothing if --without-superlu is used
    if test x"$with_superlu" != x"no" ; then
        # defaultpath
        SUPERLU_LIB_PATH="$with_superlu/$my_lib_path"
        SUPERLU_INCLUDE_PATH="$with_superlu/$my_include_path"

        # set variables so that tests can use them
        direct_SUPERLU_CPPFLAGS="-I$SUPERLU_INCLUDE_PATH -DENABLE_SUPERLU"
        SUPERLU_CPPFLAGS="-I$SUPERLU_INCLUDE_PATH -DENABLE_SUPERLU"
        CPPFLAGS="$CPPFLAGS $direct_SUPERLU_CPPFLAGS"

        # check for central header
        AC_CHECK_HEADER([$my_slu_header],
            [HAVE_SUPERLU="1"],
            [
                HAVE_SUPERLU="0"
                AC_MSG_WARN([$my_slu_header not found in $SUPERLU_INCLUDE_PATH with $CPPFLAGS])
            ])

        # if header is found check for the libs
        if test x$HAVE_SUPERLU = x1 ; then
            HAVE_SUPERLU=0

            # if neither --with-superlu-lib nor --with-superlu-blaslib was
            # given, try to link dynamically or with properly names static libs
            if test x"$with_superlu_lib$with_superlu_blaslib" = x; then
                LDFLAGS="$ac_save_LDFLAGS -L$SUPERLU_LIB_PATH"
                LIBS="$ac_save_LIBS"
                AC_CHECK_LIB([superlu], [dgssvx],
                    [
                        direct_SUPERLU_LIBS="-L$SUPERLU_LIB_PATH -lsuperlu $BLAS_LIBS $FLIBS"
                        SUPERLU_LIBS="-L$SUPERLU_LIB_PATH -lsuperlu \${BLAS_LIBS} \${FLIBS}"
                        HAVE_SUPERLU="1"
                    ], [], [$BLAS_LIBS $FLIBS])
            fi

            if test $HAVE_SUPERLU = 0 && test x"$with_superlu_lib" = x; then
                # set the default
                with_superlu_lib=superlu.a
            fi

            if test $HAVE_SUPERLU = 0 &&
                test x"$with_superlu_blaslib" = x; then
                # try system blas
                LDFLAGS="$ac_save_LDFLAGS"
                LIBS="$SUPERLU_LIB_PATH/$with_superlu_lib $BLAS_LIBS $FLIBS $ac_save_LIBS"
                AC_CHECK_FUNC([dgssvx],
                    [
                        direct_SUPERLU_LIBS="$SUPERLU_LIB_PATH/$with_superlu_lib $BLAS_LIBS $FLIBS"
                        SUPERLU_LIBS="$SUPERLU_LIB_PATH/$with_superlu_lib \${BLAS_LIBS} \${FLIBS}"
                        HAVE_SUPERLU="1"
                    ])
            fi

            # No default for with_superlu_blaslib

            if test $HAVE_SUPERLU = 0 &&
                    test x"$with_superlu_blaslib" != x; then
                # try internal blas
                LDFLAGS="$ac_save_LDFLAGS"
                LIBS="$SUPERLU_LIB_PATH/$with_superlu_lib $SUPERLU_LIB_PATH/$with_superlu_blaslib $FLIBS $ac_save_LIBS"
                AC_CHECK_FUNC([dgssvx],
                    [
                        direct_SUPERLU_LIBS="$SUPERLU_LIB_PATH/$with_superlu_lib $SUPERLU_LIB_PATH/$with_superlu_blaslib $FLIBS"
                        SUPERLU_LIBS="$SUPERLU_LIB_PATH/$with_superlu_lib $SUPERLU_LIB_PATH/$with_superlu_blaslib \${FLIBS}"
                        HAVE_SUPERLU="1"
                    ])
            fi
            
            # test whether SuperLU version is at least 4.3
            if test $HAVE_SUPERLU = "1" ; then
                AC_CHECK_DECL([SLU_DOUBLE], [SUPERLU_MIN_VERSION_4_3="1"], [], [#include <$my_slu_header>])
            fi
        fi

    else # $with_superlu = no
        HAVE_SUPERLU=0
    fi

    # Inform the user whether SuperLU was sucessfully found
    AC_MSG_CHECKING([SuperLU])
    if test x$HAVE_SUPERLU = x1 ; then
        if test "$my_slu_header" = "slu_ddefs.h"; then
            if test x$SUPERLU_MIN_VERSION_4_3 = x1 ; then
                with_superlu="yes (version 4.3 or newer)"
            else
                with_superlu="yes (version 4.2 or older, post 2005)"
            fi
        else
            with_superlu="yes (pre 2005)"
        fi
    else
        with_superlu="no"
    fi
    AC_MSG_RESULT([$with_superlu])

    # check for optional member
    if test $HAVE_SUPERLU = 1 ; then
        if test "$my_slu_header" = "slu_ddefs.h"; then
            AC_CHECK_MEMBERS([mem_usage_t.expansions],[],[],[#include "slu_ddefs.h"])
        else
            AC_CHECK_MEMBERS([mem_usage_t.expansions],[],[],[#include "dsp_defs.h"])
        fi
    fi

    # substitute variables
    if test x$HAVE_SUPERLU = x0 ; then
        SUPERLU_LIBS=
        SUPERLU_CPPFLAGS=
    fi
    AC_SUBST([SUPERLU_LIBS])
    AC_SUBST([SUPERLU_CPPFLAGS])
    #DUNE_ADD_ALL_PKG([SUPERLU], [\$(SUPERLU_CPPFLAGS)], [], [\$(SUPERLU_LIBS)])

    # tell automake
    AM_CONDITIONAL(SUPERLU, test x$HAVE_SUPERLU = x1)

    # tell the preprocessor
    if test x$HAVE_SUPERLU = x1 ; then
        AC_DEFINE([HAVE_SUPERLU], [ENABLE_SUPERLU], [Define to ENABLE_SUPERLU if SUPERLU is found])
        if test "$my_slu_header" = "slu_ddefs.h"; then
            AC_DEFINE([SUPERLU_POST_2005_VERSION], 1, [define to 1 if there is  a header slu_ddefs.h in SuperLU])
            if test x$SUPERLU_MIN_VERSION_4_3 = x1 ; then
                AC_DEFINE([SUPERLU_MIN_VERSION_4_3], 1, [define to 1 if there SLU_DOUBLE imported by header slu_ddefs.h from SuperLU])
            fi
        fi
    fi

    # summary
    #DUNE_ADD_SUMMARY_ENTRY([SuperLU],[$with_superlu])

    # restore variables
    LDFLAGS="$ac_save_LDFLAGS"
    CPPFLAGS="$ac_save_CPPFLAGS"
    LIBS="$ac_save_LIBS"
])
