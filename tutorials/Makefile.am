AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDFLAGS = $(BOOST_LDFLAGS)

LDADD =	$(top_builddir)/libopmcore.la

noinst_PROGRAMS = tutorial1
tutorial1_SOURCES   = tutorial1.cpp

if UMFPACK
noinst_PROGRAMS   += tutorial2
tutorial2_SOURCES  = tutorial2.cpp
endif

if UMFPACK
noinst_PROGRAMS   += tutorial3
tutorial3_SOURCES  = tutorial3.cpp
endif
