AM_CPPFLAGS =					\
-I$(top_srcdir)

LDADD =						\
$(top_builddir)/libopmcore.la			\
$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

noinst_PROGRAMS =				\
spu_2p						\
scaneclipsedeck

spu_2p_SOURCES = spu_2p.cpp
