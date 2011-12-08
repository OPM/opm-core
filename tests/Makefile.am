AM_CPPFLAGS =					\
-I$(top_srcdir)

LDADD =						\
$(top_builddir)/libopmcore.la			\
$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


noinst_PROGRAMS =				\
test_jacsys					\
test_readvector					\
test_sf2p

if UMFPACK
noinst_PROGRAMS += test_cfs_tpfa
endif

test_jacsys_SOURCES     = test_jacsys.cpp
test_readvector_SOURCES = test_readvector.cpp
test_sf2p_SOURCES       = test_sf2p.cpp
