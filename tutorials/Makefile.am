AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(OPM_BOOST_CPPFLAGS)

AM_LDFLAGS = $(OPM_BOOST_LDFLAGS)

LDADD =	$(top_builddir)/lib/libopmcore.la

noinst_PROGRAMS   = tutorial1
tutorial1_SOURCES = tutorial1.cpp

if UMFPACK
noinst_PROGRAMS +=				\
tutorial2					\
tutorial3					\
tutorial4

tutorial2_SOURCES = tutorial2.cpp
tutorial3_SOURCES = tutorial3.cpp
tutorial4_SOURCES = tutorial4.cpp
endif
