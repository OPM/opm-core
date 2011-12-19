AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDFLAGS = $(BOOST_LDFLAGS)

LDADD =						\
$(top_builddir)/libopmcore.la			\
$(BOOST_FILESYSTEM_LIB)				\
$(BOOST_SYSTEM_LIB)				\
$(BOOST_DATE_TIME_LIB)				\
$(BOOST_UNIT_TEST_FRAMEWORK_LIB)		\
$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

noinst_PROGRAMS =				\
scaneclipsedeck

if UMFPACK
noinst_PROGRAMS += spu_2p
spu_2p_SOURCES   = spu_2p.cpp
endif
