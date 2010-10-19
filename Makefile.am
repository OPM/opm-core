# $Id:  $

#SUBDIRS = test

ACLOCAL_AMFLAGS = -I m4

libopmpressuredir = $(includedir)/opm/pressure

lib_LTLIBRARIES = libopmpressure.la

libopmpressure_HEADERS =                        \
blas_lapack.h                                   \
coarse_conn.h                                   \
coarse_sys.h                                    \
dfs.h                                           \
flow_bc.h                                       \
fsh.h                                           \
fsh_common.h                                    \
grid.h                                          \
hash_set.h                                      \
hybsys_global.h                                 \
hybsys.h                                        \
ifsh.h                                          \
ifsh_ms.h                                       \
mimetic.h                                       \
partition.h                                     \
sparse_sys.h                                    \
well.h                                          \
GridAdapter.hpp                                 \
HybridPressureSolver.hpp


libopmpressure_la_SOURCES =			\
coarse_conn.c					\
coarse_sys.c					\
dfs.c						\
flow_bc.c					\
fsh.c						\
fsh_common.c					\
hash_set.c					\
hybsys.c					\
hybsys_global.c					\
ifsh.c						\
ifsh_ms.c					\
mimetic.c					\
partition.c					\
sparse_sys.c					\
well.c


## include $(top_srcdir)/am/global-rules
