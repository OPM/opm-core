AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDADD = $(top_builddir)/libopmcore.la 

noinst_PROGRAMS =				\
scaneclipsedeck					\
spu_2p						\
wells_example


spu_2p_SOURCES        = spu_2p.cpp
spu_2p_LDADD          =  $(LDADD) $(LIBS) $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB) $(LAPACK_LIBS) $(LIBS) $(LIBS)

wells_example_SOURCES = wells_example.cpp
