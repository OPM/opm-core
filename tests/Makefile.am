AM_CPPFLAGS =					\
-I$(top_srcdir)					\
$(BOOST_CPPFLAGS)

LDFLAGS = $(BOOST_LDFLAGS)


LDADD =	$(top_builddir)/libopmcore.la


noinst_PROGRAMS =				\
bo_fluid_p_and_z_deps                           \
bo_fluid_pressuredeps                           \
bo_fluid_test                                   \
monotcubicinterpolator_test			\
param_test					\
sparsetable_test				\
sparsevector_test				\
unit_test					\
test_readvector					\
test_sf2p

bo_fluid_p_and_z_deps_SOURCES = bo_fluid_p_and_z_deps.cpp
bo_fluid_pressuredeps_SOURCES = bo_fluid_pressuredeps.cpp
bo_fluid_test_SOURCES = bo_fluid_test.cpp
monotcubicinterpolator_test_SOURCES = monotcubicinterpolator_test.cpp
param_test_SOURCES        = param_test.cpp
sparsetable_test_SOURCES  = sparsetable_test.cpp
sparsevector_test_SOURCES = sparsevector_test.cpp
unit_test_SOURCES         = unit_test.cpp
test_readvector_SOURCES   = test_readvector.cpp
test_sf2p_SOURCES         = test_sf2p.cpp

if UMFPACK
noinst_PROGRAMS    += test_cfs_tpfa
noinst_PROGRAMS    += test_jacsys
test_jacsys_SOURCES = test_jacsys.cpp
endif
