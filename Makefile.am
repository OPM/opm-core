# $Id:  $

#SUBDIRS = test

ACLOCAL_AMFLAGS = -I m4

libopmpressuredir = $(includedir)/opm/pressure

lib_LTLIBRARIES = libopmpressure.la

libopmpressure_HEADERS =                        \
blas_lapack.h                                   \
coarse_conn.h                                   \
coarse_sys.h                                    \
compr_quant.h                                   \
dfs.h                                           \
flow_bc.h                                       \
fsh.h                                           \
grid.h                                          \
hash_set.h                                      \
hybsys_global.h                                 \
hybsys.h                                        \
ifs_tpfa.h                                      \
ifsh_ms.h                                       \
mimetic.h                                       \
partition.h                                     \
sparse_sys.h                                    \
trans_tpfa.h                                    \
well.h                                          \
GridAdapter.hpp                                 \
HybridPressureSolver.hpp


libopmpressure_la_SOURCES =                     \
cfsh.c                                          \
coarse_conn.c                                   \
coarse_sys.c                                    \
compr_quant.c                                   \
dfs.c                                           \
flow_bc.c                                       \
fsh.c                                           \
fsh_common_impl.c                               \
hash_set.c                                      \
hybsys.c                                        \
hybsys_global.c                                 \
ifs_tpfa.c                                      \
ifsh.c                                          \
ifsh_ms.c                                       \
mimetic.c                                       \
partition.c                                     \
sparse_sys.c                                    \
trans_tpfa.c                                    \
well.c


## include $(top_srcdir)/am/global-rules
