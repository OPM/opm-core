AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDFLAGS = $(BOOST_LDFLAGS)

LDADD =	$(top_builddir)/libopmcore.la

noinst_PROGRAMS =				\
scaneclipsedeck

if UMFPACK
noinst_PROGRAMS += spu_2p
spu_2p_SOURCES   = spu_2p.cpp
spu_2p_LDADD     =				\
$(BOOST_FILESYSTEM_LIB)				\
$(BOOST_SYSTEM_LIB)				\
$(LDADD)					\
$(LAPACK_LIBS)
endif
