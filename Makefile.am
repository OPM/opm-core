# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

#SUBDIRS = test

ACLOCAL_AMFLAGS = -I m4

cpgpreprocessdir = $(includedir)/dune/grid/cpgpreprocess

cpgpreprocess_HEADERS = facetopology.h preprocess.h sparsetable.h \
                     uniquepoints.h

lib_LTLIBRARIES = libcpgpreprocess.la
noinst_LTLIBRARIES = libcpgpreprocess_noinst.la

libcpgpreprocess_la_SOURCES =                   \
facetopology.c                                  \
uniquepoints.c                                  \
preprocess.c                                    \
sparsetable.c

libcpgpreprocess_noinst_la_SOURCES =            \
$(libcpgpreprocess_la_SOURCES)
