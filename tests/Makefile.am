AM_CPPFLAGS =					\
-I$(top_srcdir)

LDADD =						\
$(top_builddir)/libopmcore.la			\
$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


noinst_PROGRAMS =				\
test_readvector					\
test_sf2p

test_readvector_SOURCES = test_readvector.cpp
test_sf2p_SOURCES       = test_sf2p.cpp

if UMFPACK
noinst_PROGRAMS += test_cfs_tpfa
noinst_PROGRAMS += test_jacsys
test_jacsys_SOURCES     = test_jacsys.cpp
endif
