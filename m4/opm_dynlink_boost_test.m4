# _OPM_DYNLINK_BOOST_TEST_SRC(Symbol)
# Generate source text for use in AC_LINK_IFELSE when determining
# how to link the Boost.Test library.
#
# Note:
#   We use AC_LANG_SOURCE rather than AC_LANG_PROGRAM to avoid
#   multiple definitions of "main()" (defined by both the UTF
#   *and* the AC_LANG_PROGRAM macro).
AC_DEFUN([_OPM_DYNLINK_BOOST_TEST_SRC],
[AC_LANG_SOURCE(
[[$1
#define BOOST_TEST_MODULE OPM_DYNLINK_TEST
#include <boost/test/unit_test.hpp>

int f(int x) { return 2 * x; }

BOOST_AUTO_TEST_CASE(DynlinkConfigureTest) {
  BOOST_CHECK_MESSAGE(f(2) == 4,
                      "Apparently, multiplication doesn't "
                      "work: f(2) = " << f(2));
}
]]dnl
)[]dnl
])

dnl -------------------------------------------------------------------

# OPM_DYNLINK_BOOST_TEST
# Determine how to link (or compile+link) tests based on the UTF.
#
# If the system uses dynamic linking, then all tests need to
#
#     #define BOOST_TEST_DYN_LINK
#
# Otherwise, this symbol must *not* be #define'd.
#
# Macro defines the symbol HAVE_DYNAMIC_BOOST_TEST (to 1) if the
# system uses dynamic linking of Boost.Test .
AC_DEFUN([OPM_DYNLINK_BOOST_TEST],
[
AC_REQUIRE([AX_BOOST_BASE])
AC_REQUIRE([AX_BOOST_UNIT_TEST_FRAMEWORK])

_opm_LIBS_SAVE="$LIBS"
_opm_CPPFLAGS_SAVE="$CPPFLAGS"

LIBS="${BOOST_LDFLAGS} ${BOOST_UNIT_TEST_FRAMEWORK_LIB} ${LIBS}"
CPPFLAGS="${BOOST_CPPFLAGS} ${CPPFLAGS}"

AC_LANG_PUSH([C++])

AC_CACHE_CHECK([if the Boost.Test library can be linked statically],dnl
[opm_cv_boost_link_static],dnl
[AC_LINK_IFELSE([_OPM_DYNLINK_BOOST_TEST_SRC([])],
 [opm_cv_boost_link_static=yes],dnl
 [opm_cv_boost_link_static=no])[]dnl
])[]dnl

AC_CACHE_CHECK([if the Boost.Test library can be linked dynamically],dnl
[opm_cv_boost_link_dynamic],dnl
[AC_LINK_IFELSE([_OPM_DYNLINK_BOOST_TEST_SRC(dnl
[#define BOOST_TEST_DYN_LINK])],
 [opm_cv_boost_link_dynamic=yes],dnl
 [opm_cv_boost_link_dynamic=no])[]dnl
])[]dnl

AC_LANG_POP([C++])

LIBS="$_opm_LIBS_SAVE"
CPPFLAGS="$_opm_CPPFLAGS_SAVE"

AS_IF([test x"$opm_cv_boost_link_static" = x"yes" -o \
            x"$opm_cv_boost_link_dynamic" = x"yes"],
[AS_IF([test x"$opm_cv_boost_link_dynamic" = x"yes"],
 [AC_DEFINE([HAVE_DYNAMIC_BOOST_TEST], [1],
            [Define to `1' if Boost.Test should use BOOST_TEST_DYN_LINK])]
 [:])[]dnl
],dnl
[AC_MSG_NOTICE([Boost.Test is not supported on this system])])
])[]dnl
