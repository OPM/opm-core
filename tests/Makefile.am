AM_CPPFLAGS =					\
-I$(top_srcdir)                                 \
$(BOOST_CPPFLAGS)

LDFLAGS = $(BOOST_LDADD)


LDADD =						\
$(top_builddir)/libopmcore.la			\
$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)    \
$(BOOST_FILESYSTEM_LIB)                         \
$(BOOST_SYSTEM_LIB)                             \
$(BOOST_DATE_TIME_LIB)                          \
$(BOOST_UNIT_TEST_FRAMEWORK_LIB)


noinst_PROGRAMS =				\
monotcubicinterpolator_test                     \
sparsetable_test                                \
sparsevector_test                               \
unit_test                                       \
test_readvector					\
test_sf2p

monotcubicinterpolator_test_SOURCES = monotcubicinterpolator_test.cpp
sparsetable_test_SOURCES = sparsetable_test.cpp
sparsevector_test_SOURCES = sparsevector_test.cpp
unit_test_SOURCES       = unit_test.cpp
test_readvector_SOURCES = test_readvector.cpp
test_sf2p_SOURCES       = test_sf2p.cpp

if UMFPACK
noinst_PROGRAMS += test_cfs_tpfa
noinst_PROGRAMS += test_jacsys
test_jacsys_SOURCES     = test_jacsys.cpp
endif
