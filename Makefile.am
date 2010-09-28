# $Id:  $

#SUBDIRS = test

ACLOCAL_AMFLAGS = -I m4

libmimeticdir = $(includedir)/mimetic

lib_LTLIBRARIES = libmimetic.la

libmimetic_HEADERS =                            \
blas_lapack.h                                   \
coarse_conn.h                                   \
coarse_sys.h                                    \
dfs.h                                           \
flow_bc.h                                       \
grid.h                                          \
hash_set.h                                      \
hybsys_global.h                                 \
hybsys.h                                        \
ifsh.h                                          \
mimetic.h                                       \
partition.h                                     \
sparse_sys.h                                    \
well.h


libmimetic_la_SOURCES =                         \
coarse_conn.c                                   \
coarse_sys.c                                    \
dfs.c                                           \
flow_bc.c                                       \
hash_set.c                                      \
hybsys.c                                        \
hybsys_global.c                                 \
ifsh.c                                          \
mimetic.c                                       \
partition.c                                     \
sparse_sys.c                                    \
well.c


## include $(top_srcdir)/am/global-rules
