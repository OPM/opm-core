# Additional definitions to ensure proper include paths

AM_CPPFLAGS =					\
-I$(top_srcdir)/src

AM_LDFLAGS =					\
$(LAPACK_LIBS) $(BLAS_LIBS)

#-----------------------------------------------------------------------

# Declare example programs
noinst_PROGRAMS =

# "spu_2p" depends on availability of OPMPressure and UMFPACK
if UMFPACK
if OPMPRESSURE_SPARSE_SYS
noinst_PROGRAMS += spu_2p
endif
endif

#-----------------------------------------------------------------------

# Example program sources

spu_2p_SOURCES =				\
spu_2p.cpp					\
call_umfpack.c					\
cart_grid.c					\
transport_source.c
