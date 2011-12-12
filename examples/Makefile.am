AM_CPPFLAGS =					\
-I$(top_srcdir)

LDADD =						\
$(top_builddir)/libopmcore.la			\
$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

noinst_PROGRAMS =				\
scaneclipsedeck

if UMFPACK
noinst_PROGRAMS += spu_2p
spu_2p_SOURCES = spu_2p.cpp
endif
