# Build-time flags needed to form example programs
AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

# All targets link to the library
LDADD =						\
$(top_builddir)/lib/libopmcore.la		\
$(BOOST_FILESYSTEM_LIB)				\
$(BOOST_SYSTEM_LIB)

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

refine_wells_SOURCES          = refine_wells.cpp
sim_2p_incomp_reorder_SOURCES = sim_2p_incomp_reorder.cpp
sim_wateroil_SOURCES          = sim_wateroil.cpp
wells_example_SOURCES         = wells_example.cpp

# ----------------------------------------------------------------------
# Optional examples, or examples that use optional add-on components.

if UMFPACK
noinst_PROGRAMS += spu_2p

spu_2p_SOURCES   = spu_2p.cpp
spu_2p_LDADD     =				\
$(LDADD)					\
$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
endif
