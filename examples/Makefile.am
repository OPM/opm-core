AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDADD = $(top_builddir)/libopmcore.la

noinst_PROGRAMS =				\
refine_wells					\
scaneclipsedeck					\
sim_wateroil					\
wells_example					\
sim_2p_incomp_reorder

refine_wells_SOURCES  = refine_wells.cpp

sim_wateroil_SOURCES  = sim_wateroil.cpp
sim_wateroil_LDADD    =				\
$(LDADD) $(LIBS)				\
$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)	\
$(LAPACK_LIBS) $(LIBS) $(LIBS)

sim_2p_incomp_reorder_SOURCES  = sim_2p_incomp_reorder.cpp
sim_2p_incomp_reorder_LDADD    =		\
$(LDADD) $(LIBS)				\
$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)	\
$(LAPACK_LIBS) $(LIBS) $(LIBS)

wells_example_SOURCES = wells_example.cpp

if UMFPACK
noinst_PROGRAMS += spu_2p

spu_2p_SOURCES = spu_2p.cpp
spu_2p_LDADD   =				\
$(LDADD) $(LIBS)				\
$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)	\
$(LAPACK_LIBS) $(LIBS) $(LIBS)
endif
