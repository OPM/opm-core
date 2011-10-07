# $Date$
# $Revision$

check_PROGRAMS = sparsetable_test sparsevector_test monotcubicinterpolator_test
noinst_PROGRAMS = unit_test

sparsetable_test_SOURCES = sparsetable_test.cpp
sparsetable_test_CXXFLAGS = $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS)
sparsetable_test_LDADD = $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) \
                         $(DUNE_LIBS) $(DUNEMPILIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

sparsevector_test_SOURCES = sparsevector_test.cpp
sparsevector_test_CXXFLAGS = $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS)
sparsevector_test_LDADD = $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) \
                          $(DUNE_LIBS) $(DUNEMPILIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

monotcubicinterpolator_test_SOURCES = monotcubicinterpolator_test.cpp
monotcubicinterpolator_test_CXXFLAGS = $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS)
monotcubicinterpolator_test_LDADD = $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS) \
                                    $(DUNE_LIBS) $(DUNEMPILIBS) ../libcommon.la

unit_test_SOURCES = unit_test.cpp
unit_test_CXXFLAGS = $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS)
unit_test_LDADD = $(DUNE_LDFLAGS) $(DUNEMPILDFLAGS) $(DUNE_LIBS) $(DUNEMPILIBS)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules

