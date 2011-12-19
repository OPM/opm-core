# $Date$
# $Revision$

SUBDIRS = tinyxml . test

paramdir = $(includedir)/opm/core/utility/parameters/param

param_HEADERS = Parameter.hpp ParameterGroup.hpp ParameterGroup_impl.hpp \
                ParameterMapItem.hpp ParameterRequirement.hpp \
                ParameterStrings.hpp ParameterTools.hpp ParameterXML.hpp

noinst_LTLIBRARIES = libparam.la

libparam_la_SOURCES = Parameter.cpp ParameterGroup.cpp ParameterTools.cpp \
                      ParameterXML.cpp

libparam_la_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(BOOST_CPPFLAGS)
libparam_la_LDFLAGS  = $(ALL_PKG_LDFLAGS) $(BOOST_LDFLAGS)
libparam_la_LIBADD = tinyxml/libtinyxml.la  $(ALL_PKG_LIBS) $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)

include $(top_srcdir)/am/global-rules
