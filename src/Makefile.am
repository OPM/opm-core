libopmpressuredir = $(includedir)/opm/pressure

lib_LTLIBRARIES    = libopmpressure.la
noinst_LTLIBRARIES = libopmpressure_noinst.la

libopmpressure_HEADERS =                        \
blas_lapack.h                                   \
coarse_conn.h                                   \
coarse_sys.h                                    \
compr_quant.h                                   \
compr_quant_general.h                           \
cfs_tpfa.h                                      \
cfs_tpfa_residual.h                             \
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
HybridPressureSolver.hpp                        \
TPFAPressureSolver.hpp


libopmpressure_la_SOURCES =                     \
cfsh.c                                          \
coarse_conn.c                                   \
coarse_sys.c                                    \
compr_quant.c                                   \
compr_quant_general.c                           \
cfs_tpfa.c                                      \
cfs_tpfa_residual.c                             \
dfs.c                                           \
flow_bc.c                                       \
fsh.c                                           \
fsh_common_impl.h                               \
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

libopmpressure_noinst_la_SOURCES = $(libopmpressure_la_SOURCES)


## include $(top_srcdir)/am/global-rules
