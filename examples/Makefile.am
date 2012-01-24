AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDFLAGS = $(BOOST_LDFLAGS)

LDADD =	$(top_builddir)/libopmcore.la

noinst_PROGRAMS =				\
scaneclipsedeck

if UMFPACK
noinst_PROGRAMS += spu_2p reorder-qfs

spu_2p_SOURCES   = spu_2p.cpp
reorder_qfs_SOURCES = reorder-qfs.cpp
endif
