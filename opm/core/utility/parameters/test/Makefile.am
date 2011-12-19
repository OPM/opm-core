# $Date: 2009-12-08 15:21:36 +0100 (Tue, 08 Dec 2009) $
# $Revision: 699 $

check_PROGRAMS = param_test

param_test_SOURCES = param_test.cpp
param_test_CXXFLAGS = $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS)
param_test_LDADD = $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS) $(DUNEMPILIBS) $(BOOST_LDFLAGS) \
                         $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) ../libparam.la

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

