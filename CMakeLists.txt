# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:
cmake_minimum_required (VERSION 2.8)
project (opm-core)
set (opm-core_VERSION_MAJOR 0)
set (opm-core_VERSION_MINOR 3)
enable_language (CXX)

# build debug by default
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Debug")
endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# all public header files are together with the source
set (opm-core_INCLUDE_DIR "${PROJECT_SOURCE_DIR}")
list (APPEND opm-core_INCLUDE_DIRS "${opm-core_INCLUDE_DIR}")

# additional search modules
set (opm-core_MODULE_DIR "${PROJECT_SOURCE_DIR}/cmake/Modules")
list (APPEND CMAKE_MODULE_PATH ${opm-core_MODULE_DIR})

# use tricks to do faster builds
include (UseFastBuilds)

# macro to set standard variables (INCLUDE_DIRS, LIBRARIES etc.)
include (OpmFind)

# dependencies
list (APPEND opm-core_DEPS
	# compile with C++0x/11 support if available
	"CXX11Features"
	# matrix library
	"BLAS REQUIRED"
	"LAPACK REQUIRED"
	# Tim Davis' SuiteSparse archive	
	"SuiteSparse COMPONENTS umfpack"
	# solver
	"SUPERLU"
	# xml processing (for config parsing)	
	"LibXml2 REQUIRED"
	# various runtime library enhancements
	"Boost 1.39.0
		COMPONENTS date_time filesystem system unit_test_framework REQUIRED"
	# DUNE dependency
	"dune-istl"
	)
find_and_append_package_list (${opm-core_DEPS})

# put debug information into every executable
include (UseDebugSymbols)

# optimize full if we're not doing a debug build
include (UseOptimization)

# turn on all warnings
include (UseWarnings)

# detect if Boost is in a shared library
include (UseDynamicBoost)

# put libraries in lib/
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")

# find all the source code (note that these variables have name after
# the target library and not the project)
file (GLOB_RECURSE opmcore_SOURCES "opm/core/*.c" "opm/core/*.cpp")
file (GLOB_RECURSE opmcore_HEADERS "opm/core/*.h" "opm/core/*.hpp")

# these files are provided in source control, but can only compile with Matlab
# available
list (REMOVE_ITEM opmcore_SOURCES
	${PROJECT_SOURCE_DIR}/opm/core/grid/cpgpreprocess/mxgrdecl.c
	${PROJECT_SOURCE_DIR}/opm/core/grid/cpgpreprocess/processgrid.c
	${PROJECT_SOURCE_DIR}/opm/core/utility/parameters/tinyxml/xmltest.cpp
	${PROJECT_SOURCE_DIR}/opm/core/linalg/LinearSolverAGMG.cpp
	)

# we don't try to find any of these
set (HAVE_AGMG 0)
set (HAVE_ERT)

# create configuration header which describes available features
# necessary to compile this library. singular version is the names that
# is required by this project alone, plural version transitively
# includes the necessary defines by the dependencies
include (ConfigVars)
set (opm-core_CONFIG_VAR
	HAVE_AGMG
	HAVE_DUNE_ISTL
	HAVE_DYNAMIC_BOOST_TEST
	HAVE_ERT
	HAVE_SUITESPARSE_UMFPACK_H
	HAVE_NULLPTR
	HAVE_STATIC_ASSERT
	)
list (APPEND opm-core_CONFIG_VARS ${opm-core_CONFIG_VAR})
configure_vars (
	FILE  CXX  "${PROJECT_BINARY_DIR}/config.h"
	WRITE ${opm-core_CONFIG_VARS}
	)

# some CMake properties do not do list expansion
list (REMOVE_DUPLICATES opm-core_LINKER_FLAGS)
string (REPLACE ";" " " opm-core_LINKER_FLAGS_STR "${opm-core_LINKER_FLAGS}")
remove_duplicate_libraries (opm-core)

# create this library
include_directories (${opm-core_INCLUDE_DIRS})
link_directories (${opm-core_LIBRARY_DIRS})
if (opm-core_DEFINITIONS)
	list (REMOVE_DUPLICATES opm-core_DEFINITIONS)
	add_definitions (${opm-core_DEFINITIONS})
endif (opm-core_DEFINITIONS)
add_library (opmcore SHARED ${opmcore_SOURCES})
set (opm-core_VERSION "${opm-core_VERSION_MAJOR}.${opm-core_VERSION_MINOR}")
set_target_properties (opmcore PROPERTIES
	SOVERSION ${opm-core_VERSION_MAJOR}
	VERSION ${opm-core_VERSION}
	LINK_FLAGS "${opm-core_LINKER_FLAGS_STR}"
	)
target_link_libraries (opmcore ${opm-core_LIBRARIES})

# queue this executable to be stripped
strip_debug_symbols (LIBRARY opmcore)

# we need to know the name of the library which is generated
get_target_property (opm-core_LIBRARY opmcore LOCATION)

# write configuration file to locate library
configure_file (
	${PROJECT_SOURCE_DIR}/opm-core-config.cmake.in
	${PROJECT_BINARY_DIR}/opm-core-config.cmake
	@ONLY
	)
configure_vars (
	FILE CMAKE "${PROJECT_BINARY_DIR}/opm-core-config.cmake"
	APPEND "${opm-core_CONFIG_VARS}"
	)
include (WriteBasicConfigVersionFile)
write_basic_config_version_file (
	"${PROJECT_BINARY_DIR}/opm-core-config-version.cmake"
	VERSION ${opm-core_VERSION}
	COMPATIBILITY SameMajorVersion
	)

### installation ###
string (LENGTH ${PROJECT_SOURCE_DIR} _prefix_length)
foreach (_hdr IN LISTS opmcore_HEADERS)
	get_filename_component (_dir ${_hdr} PATH)
	string (SUBSTRING ${_dir} ${_prefix_length} -1 _rel_dir)
	install (
		FILES ${_hdr}
		DESTINATION include${_rel_dir}
		)
endforeach (_hdr)
install (
	TARGETS opmcore
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)
message (STATUS "This build defaults to installing in ${CMAKE_INSTALL_PREFIX}")

# installation of CMake modules to help user programs locate the library
include (OpmProject)
opm_cmake_config (opm-core)

### test programs ###

# find the source code
file (GLOB_RECURSE tests_SOURCES "tests/*.c" "tests/*.cpp")

# conditionally disable tests when features aren't available
macro (cond_disable_test name)
	if ((NOT DEFINED ${HAVE_${name}}) OR (NOT ${HAVE_${name}}))
		message (STATUS "${name} test disabled, since ${name} is not found.")
		string (TOLOWER "${name}" name_lower)
		get_filename_component (test_${name}_FILE "tests/test_${name_lower}.cpp" ABSOLUTE)
		list (REMOVE_ITEM tests_SOURCES "${test_${name}_FILE}")
	endif ((NOT DEFINED ${HAVE_${name}}) OR (NOT ${HAVE_${name}}))
endmacro (cond_disable_test name)
cond_disable_test ("AGMG")
cond_disable_test ("ERT")

# compile each of these separately
foreach (test_FILE IN LISTS tests_SOURCES)
	get_filename_component (test_NAME "${test_FILE}" NAME_WE)
	add_executable (${test_NAME} ${test_FILE})
	set_target_properties (${test_NAME} PROPERTIES
		LINK_FLAGS "${opm-core_LINKER_FLAGS_STR}"
		)
	target_link_libraries (${test_NAME} opmcore ${opm-core_LIBRARIES})
	strip_debug_symbols (RUNTIME ${test_NAME})
endforeach (test_FILE)
