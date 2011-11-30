# $Id: duneproject 5489 2009-03-25 11:19:24Z sander $

#SUBDIRS = test

ACLOCAL_AMFLAGS = -I m4

cpgpreprocessdir = $(includedir)/cpgpreprocess

cpgpreprocess_HEADERS = cgridinterface.h preprocess.h


lib_LTLIBRARIES = libcpgpreprocess.la
noinst_LTLIBRARIES = libcpgpreprocess_noinst.la

libcpgpreprocess_la_SOURCES =			\
cgridinterface.c				\
facetopology.c					\
uniquepoints.c					\
preprocess.c					\
sparsetable.c					\
facetopology.h					\
uniquepoints.h					\
sparsetable.h

libcpgpreprocess_noinst_la_SOURCES =            \
$(libcpgpreprocess_la_SOURCES)
