# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:
cmake_minimum_required (VERSION 2.8)
project (opm-core)
set (opm-core_MAJOR_VERSION 0)
set (opm-core_MINOR_VERSION 3)
enable_language (CXX)

# all public header files are together with the source
set (opm-core_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}")

# additional search modules
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

# put libraries in lib/
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

# find all the source code
file (GLOB_RECURSE opm-core_SOURCES "opm/*.c" "opm/*.cpp")

# these files are provided in source control, but can only compile with Matlab
# available
list (REMOVE_ITEM opm-core_SOURCES
	${PROJECT_SOURCE_DIR}/opm/core/grid/cpgpreprocess/mxgrdecl.c
	${PROJECT_SOURCE_DIR}/opm/core/grid/cpgpreprocess/processgrid.c
	${PROJECT_SOURCE_DIR}/opm/core/utility/parameters/tinyxml/xmltest.cpp
	${PROJECT_SOURCE_DIR}/opm/core/linalg/LinearSolverAGMG.cpp
	)

# we don't try to find any of these
set (HAVE_AGMG 0)
set (HAVE_DUNE_ISTL 0)
set (HAVE_DYNAMIC_BOOST_TEST 0)
set (HAVE_ERT)
set (HAVE_SUITESPARSE_UMFPACK_H 0)

# create configuration header which describes available features
# necessary to compile this library
include (UseConfigVars)
list (APPEND opm-core_CONFIG_VARS
	"HAVE_AGMG"
	"HAVE_DUNE_ISTL"
	"HAVE_DYNAMIC_BOOST_TEST"
	"HAVE_ERT"
	"HAVE_SUITESPARSE_UMFPACK_H"
	)
configure_vars (
	FILE  "${PROJECT_BINARY_DIR}/config.h"
	WRITE ${opm-core_CONFIG_VARS}
	)

# create this library
include_directories (${opm-core_INCLUDE_DIRS})
add_library (opmcore SHARED ${opm-core_SOURCES})
set_target_properties (opmcore PROPERTIES
	SOVERSION ${opm-core_MAJOR_VERSION}
	VERSION ${opm-core_MAJOR_VERSION}.${opm-core_MINOR_VERSION}
	)
