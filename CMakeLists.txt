# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# key information about the library
set (opm-core_NAME "opm-core")
set (opm-core_DESCRIPTION "Open Porous Media Initiative Core Library")
set (opm-core_DIR "opm/core")
set (opm-core_VERSION_MAJOR 1)
set (opm-core_VERSION_MINOR 0)

# C++ project
cmake_minimum_required (VERSION 2.8)
project (${opm-core_NAME})
enable_language (C)
enable_language (CXX)

# additional search modules
set (opm-core_MODULE_DIR "${PROJECT_SOURCE_DIR}/cmake/Modules")
list (APPEND CMAKE_MODULE_PATH ${opm-core_MODULE_DIR})

# print system information to better pinpoint issues from log alone
include (UseSystemInfo)
system_info ()

# very early try to print repo id (to pinpoint version if something goes wrong)
include (UseVCSInfo)
vcs_info ()

# include special
if (CMAKE_VERSION VERSION_LESS "2.8.7")
	message (STATUS "Enabling backward compatibility modules for CMake ${CMAKE_VERSION}")
	list (APPEND CMAKE_MODULE_PATH "${opm-core_MODULE_DIR}/compat-2.8.7")
endif (CMAKE_VERSION VERSION_LESS "2.8.7")

# default settings: build static debug library
include (OpmDefaults)
opm_defaults (opm-core)
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# use tricks to do faster builds
include (UseFastBuilds)

# precompiled headers
include (UsePrecompHeaders)

# macro to set standard variables (INCLUDE_DIRS, LIBRARIES etc.)
include (OpmFind)

# dependencies
list (APPEND opm-core_DEPS
	# compile with C99 support if available
	"C99"
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
	"TinyXML"
	# various runtime library enhancements
	"Boost 1.39.0
		COMPONENTS date_time filesystem system unit_test_framework REQUIRED"
	# DUNE dependency
	"dune-istl"
	# Ensembles-based Reservoir Tools (ERT)
	"ERT"
	)
find_and_append_package_list_to (opm-core ${opm-core_DEPS})

# remove the dependency on the testing framework from the main library;
# it is not possible to query for Boost twice with different components.
list (REMOVE_ITEM opm-core_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# don't import more libraries than we need to
include (UseOnlyNeeded)

# put debug information into every executable
include (UseDebugSymbols)

# optimize full if we're not doing a debug build
include (UseOptimization)

# turn on all warnings
include (UseWarnings)

# detect if Boost is in a shared library
include (UseDynamicBoost)

# needed for Debian installation scheme
include (UseMultiArch)

# this module contains code to figure out which files is where
include (OpmFiles)
opm_auto_dirs ()

# put libraries in lib/
opm_out_dirs ()

# identify the compilation units in the library
opm_sources (opm-core)

# enumerate all testing programs in test/ directory
opm_find_tests ()

# tutorial programs are found in the tutorials/ directory
opm_find_tutorials ()

# example programs are found in the examples/ directory
opm_find_examples ()

# Algebraic Multigrid must be compiled together with our program;
# if it is not available, then remove our corresponding component
find_package (AGMG)
if (AGMG_FOUND)
	list (APPEND opm-core_SOURCES ${AGMG_SOURCES})
endif (AGMG_FOUND)

# these solvers are only compiled in if their dependency is found
if (NOT AGMG_FOUND)
	list (REMOVE_ITEM opm-core_SOURCES
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/linalg/LinearSolverAGMG.cpp
		)
endif (NOT AGMG_FOUND)
if (NOT dune-istl_FOUND)
	list (REMOVE_ITEM opm-core_SOURCES
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/linalg/LinearSolverIstl.cpp
		)
endif (NOT dune-istl_FOUND)
if (NOT SuiteSparse_FOUND)
	list (REMOVE_ITEM opm-core_SOURCES
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/linalg/call_umfpack.c
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/linalg/LinearSolverUmfpack.cpp
		)
	list (REMOVE_ITEM tutorial_SOURCES
		${PROJECT_SOURCE_DIR}/${tutorial_DIR}/tutorial2.cpp
		${PROJECT_SOURCE_DIR}/${tutorial_DIR}/tutorial3.cpp
		${PROJECT_SOURCE_DIR}/${tutorial_DIR}/tutorial4.cpp
		)
	list (REMOVE_ITEM examples_SOURCES
		${PROJECT_SOURCE_DIR}/${examples_DIR}/spu_2p.cpp
		)
endif (NOT SuiteSparse_FOUND)

# these files are provided in source control, but can only compile with Matlab
# available; we are not supposed to include the TinyXML test prog. regardless
list (REMOVE_ITEM opm-core_SOURCES
	${PROJECT_SOURCE_DIR}/${opm-core_DIR}/grid/cpgpreprocess/mxgrdecl.c
	${PROJECT_SOURCE_DIR}/${opm-core_DIR}/grid/cpgpreprocess/processgrid.c
	${PROJECT_SOURCE_DIR}/${opm-core_DIR}/utility/parameters/tinyxml/xmltest.cpp
	)

# remove inline TinyXML if a system version was found
if (TinyXML_FOUND)
	file (GLOB_RECURSE _inline_tinyxml "${opm-core_DIR}/utility/parameters/tinyxml/*")
	foreach (_file IN LISTS _inline_tinyxml)
		list (REMOVE_ITEM opm-core_SOURCES ${_file})
	endforeach (_file)
endif (TinyXML_FOUND)

# anyhow remove it from the header list (so it doesn't get installed)
list (REMOVE_ITEM opm-core_HEADERS "${opm-core_DIR}/utility/parameters/tinyxml/tinystr.h")
list (REMOVE_ITEM opm-core_HEADERS "${opm-core_DIR}/utility/parameters/tinyxml/tinyxml.h")

# HAVE_ERT is used as an #ifdef, not as an #if in the source code, if it
# is not true, then it should be unset altogether
if (NOT HAVE_ERT)
	set (HAVE_ERT)
	list (REMOVE_ITEM examples_SOURCES
		${PROJECT_SOURCE_DIR}/examples/import_rewrite.cpp
		)
endif (NOT HAVE_ERT)

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
include (UseFortranWrappers)
define_fc_func (
	APPEND "${PROJECT_BINARY_DIR}/config.h"
	)

# compile main library; pull in all required includes and libraries
include (OpmCompile)
opm_compile (opm-core)

# installation target: copy the library together with debug and
# configuration files to system directories
include (OpmInstall)
opm_install (opm-core)
message (STATUS "This build defaults to installing in ${CMAKE_INSTALL_PREFIX}")

# installation of CMake modules to help user programs locate the library
include (OpmProject)
opm_cmake_config (opm-core)

# routines to build satellites such as tests, tutorials and samples
include (OpmSatellites)

# tutorial programs are found in the tutorials/ directory
opm_compile_satellites (opm-core tutorial "" "")
opm_compile_satellites (opm-core examples "" "")

# infrastructure for testing
enable_testing ()
include (CTest)

# conditionally disable tests when features aren't available
macro (cond_disable_test name)
	if ((NOT DEFINED HAVE_${name}) OR (NOT HAVE_${name}))
		message (STATUS "${name} test disabled, since ${name} is not found.")
		string (TOLOWER "${name}" name_lower)
		get_filename_component (test_${name}_FILE "tests/test_${name_lower}.cpp" ABSOLUTE)
		list (REMOVE_ITEM tests_SOURCES "${test_${name}_FILE}")
	endif ((NOT DEFINED HAVE_${name}) OR (NOT HAVE_${name}))
endmacro (cond_disable_test name)
cond_disable_test ("AGMG")
cond_disable_test ("ERT")

# make datafiles necessary for tests available in output directory
opm_data (tests datafiles "tests/*.xml")
opm_compile_satellites (opm-core tests "" "^test_([^/]*)$")

# use this target to run all tests
add_custom_target (check
	COMMAND ${CMAKE_CTEST_COMMAND}
	DEPENDS tests
	COMMENT "Checking if library is functional"
	VERBATIM
	)

# generate documentation from source code with Doxygen;
# setup install target for this documentation
set (doxy_dir "Documentation")
include (OpmDoc)
opm_doc (opm-core ${doxy_dir})

# provide compatibility with using this build in dunecontrol
include (DuneCompat)
include (LibtoolArchives)
configure_la (opm-core ${opm-core_TARGET} opm-core_LIBTOOL_ARCHIVE)

### clean in-source builds ###
include (OpmDistClean)
opm_dist_clean (opm-core)

# smart wrapper that auto-parallelizes builds
file (COPY
	GNUmakefile
	DESTINATION ${PROJECT_BINARY_DIR}
	)
