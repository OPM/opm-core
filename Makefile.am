# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

#SUBDIRS = test

preprocessdir = $(includedir)/dune/grid/preprocess

preprocess_HEADERS = facetopology.h preprocess.h sparsetable.h \
                     uniquepoints.h

lib_LTLIBRARIES = libpreprocess.la

libpreprocess_la_SOURCES = facetopology.c  uniquepoints.c \
			   matalloc.c      preprocess.c   \
                           sparsetable.c


include $(top_srcdir)/am/global-rules

