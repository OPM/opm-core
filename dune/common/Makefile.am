# $Date$
# $Revision$

SUBDIRS = . param test

common_HEADERS	= \
	Average.hpp \
	EclipseGridInspector.hpp \
	EclipseGridParser.hpp \
	EclipseGridParserHelpers.hpp \
        ErrorMacros.hpp \
	Factory.hpp \
	linInt.hpp \
        MonotCubicInterpolator.hpp \
	RootFinders.hpp \
	SparseTable.hpp \
	SparseVector.hpp \
	SpecialEclipseFields.hpp \
	StopWatch.hpp \
	Units.hpp

commondir = $(includedir)/dune/common

noinst_LTLIBRARIES = libcommon.la

libcommon_la_SOURCES = \
	EclipseGridInspector.cpp \
	EclipseGridParser.cpp \
	MonotCubicInterpolator.cpp \
	StopWatch.cpp


libcommon_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(BOOST_CPPFLAGS)
libcommon_la_LDFLAGS  = $(DUNE_LDFLAGS) $(BOOST_LDFLAGS) \
                        $(BOOST_DATE_TIME_LIB)


include $(top_srcdir)/am/global-rules
