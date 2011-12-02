# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

#SUBDIRS = test

ACLOCAL_AMFLAGS = -I m4

cpgpreprocessdir = $(includedir)/cpgpreprocess

cpgpreprocess_HEADERS = cgridinterface.h preprocess.h

noinst_PROGRAMS = grid_io_demo

lib_LTLIBRARIES = libcpgpreprocess.la
noinst_LTLIBRARIES = libcpgpreprocess_noinst.la

libcpgpreprocess_la_SOURCES =			\
cgridinterface.c				\
facetopology.c					\
facetopology.h					\
geometry.c					\
geometry.h					\
grid.h						\
preprocess.c					\
readvector.cpp					\
readvector.hpp					\
sparsetable.c					\
sparsetable.h					\
uniquepoints.c					\
uniquepoints.h

libcpgpreprocess_noinst_la_SOURCES =            \
$(libcpgpreprocess_la_SOURCES)

grid_io_demo_SOURCES = test_readvector.cpp
grid_io_demo_LDADD = libcpgpreprocess.la