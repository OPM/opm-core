# Build-time flags needed to form example programs
AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

# All targets link to the library
LDADD = $(top_builddir)/libopmcore.la

# Convenience make variable to avoid having to list the entire
# link-time dependency for every target that requires more than
# libopmcore.la
FULL_LDADD    =					\
$(LDADD) $(LIBS)				\
$(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)	\
$(LAPACK_LIBS) $(LIBS) $(LIBS)

# ----------------------------------------------------------------------
# Declare products (i.e., the example programs).
#
# Please keep the list sorted.

noinst_PROGRAMS =				\
refine_wells					\
scaneclipsedeck					\
sim_2p_incomp_reorder				\
sim_wateroil					\
wells_example

# ----------------------------------------------------------------------
# Product constituents.  Must be specified for every product that's
# built from more than a single ".c" file and/or that link to anything
# more than the OPM-Core library.
#
# Please maintain sort order from "noinst_PROGRAMS".

refine_wells_SOURCES  = refine_wells.cpp

sim_2p_incomp_reorder_SOURCES  = sim_2p_incomp_reorder.cpp
sim_2p_incomp_reorder_LDADD    = $(FULL_LDADD)

sim_wateroil_SOURCES  = sim_wateroil.cpp
sim_wateroil_LDADD    = $(FULL_LDADD)

wells_example_SOURCES = wells_example.cpp

# ----------------------------------------------------------------------
# Optional examples, or examples that use optional add-on components.

if UMFPACK
noinst_PROGRAMS += spu_2p

spu_2p_SOURCES = spu_2p.cpp
spu_2p_LDADD   = $(FULL_LDADD)
endif
