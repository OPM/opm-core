AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDADD = $(top_builddir)/libopmcore.la

FULL_LDADD    =					\
$(LDADD) $(LIBS)				\
$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)	\
$(LAPACK_LIBS) $(LIBS) $(LIBS)



noinst_PROGRAMS =				\
refine_wells					\
scaneclipsedeck					\
sim_2p_incomp_reorder				\
sim_wateroil					\
wells_example

refine_wells_SOURCES  = refine_wells.cpp

sim_2p_incomp_reorder_SOURCES  = sim_2p_incomp_reorder.cpp
sim_2p_incomp_reorder_LDADD    = $(FULL_LDADD)

sim_wateroil_SOURCES  = sim_wateroil.cpp
sim_wateroil_LDADD    = $(FULL_LDADD)

wells_example_SOURCES = wells_example.cpp

if UMFPACK
noinst_PROGRAMS += spu_2p

spu_2p_SOURCES = spu_2p.cpp
spu_2p_LDADD   = $(FULL_LDADD)
endif
