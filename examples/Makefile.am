AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDADD =	$(top_builddir)/libopmcore.la

noinst_PROGRAMS =				\
scaneclipsedeck

if UMFPACK
noinst_PROGRAMS +=				\
spu_2p						\
wells_example

spu_2p_SOURCES        = spu_2p.cpp
wells_example_SOURCES = wells_example.cpp
endif
