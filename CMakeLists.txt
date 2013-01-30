# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:
cmake_minimum_required (VERSION 2.8)
project (opm-core)
set (opm-core_VERSION_MAJOR 1)
set (opm-core_VERSION_MINOR 0)
enable_language (C)
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

# include special
if (CMAKE_VERSION VERSION_LESS "2.8.7")
	message (STATUS "Enabling backward compatibility modules for CMake ${CMAKE_VERSION}")
	list (APPEND CMAKE_MODULE_PATH "${opm-core_MODULE_DIR}/compat-2.8.7")
endif (CMAKE_VERSION VERSION_LESS "2.8.7")

# don't import more libraries than we need to
include (UseOnlyNeeded)

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
find_and_append_package_list (${opm-core_DEPS})

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

# put libraries in lib/ (no multi-arch support in build tree)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set (CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/CMakeFiles")

# find all the source code (note that these variables have name after
# the target library and not the project). the documentation recommends
# against using globs to enumerate source code, but if you list the
# files explicitly you'll change the build files every time you add to
# the project as well as having to rebuild completely anyway.
file (GLOB_RECURSE opmcore_C_SOURCES "opm/core/[^.]*.c")
file (GLOB_RECURSE opmcore_CXX_SOURCES "opm/core/[^.]*.cpp")
file (GLOB_RECURSE opmcore_C_HEADERS "opm/core/[^.]*.h")
file (GLOB_RECURSE opmcore_CXX_HEADERS "opm/core/[^.]*.hpp")

# remove pre-compile headers from output list
set (opmcore_PRECOMP_CXX_HEADER "opm/core/opm-core-pch.hpp")
list (REMOVE_ITEM opmcore_CXX_HEADERS
	${PROJECT_SOURCE_DIR}/${opmcore_PRECOMP_CXX_HEADER}
	)

# merge both languages into one compilation/installation
set (opmcore_SOURCES ${opmcore_C_SOURCES} ${opmcore_CXX_SOURCES})
set (opmcore_HEADERS ${opmcore_C_HEADERS} ${opmcore_CXX_HEADERS})

# Algebraic Multigrid must be compiled together with our program;
# if it is not available, then remove our corresponding component
find_package (AGMG)
if (AGMG_FOUND)
	list (APPEND opmcore_SOURCES ${AGMG_SOURCES})
endif (AGMG_FOUND)

# these solvers are only compiled in if their dependency is found
if (NOT AGMG_FOUND)
	list (REMOVE_ITEM opmcore_SOURCES
		${PROJECT_SOURCE_DIR}/opm/core/linalg/LinearSolverAGMG.cpp
		)
endif (NOT AGMG_FOUND)
if (NOT dune-istl_FOUND)
	list (REMOVE_ITEM opmcore_SOURCES
		${PROJECT_SOURCE_DIR}/opm/core/linalg/LinearSolverIstl.cpp
		)
endif (NOT dune-istl_FOUND)
if (NOT SuiteSparse_FOUND)
	list (REMOVE_ITEM opmcore_SOURCES
		${PROJECT_SOURCE_DIR}/opm/core/linalg/call_umfpack.c
		${PROJECT_SOURCE_DIR}/opm/core/linalg/LinearSolverUmfpack.cpp
		)
endif (NOT SuiteSparse_FOUND)

# these files are provided in source control, but can only compile with Matlab
# available; we are not supposed to include the TinyXML test prog. regardless
list (REMOVE_ITEM opmcore_SOURCES
	${PROJECT_SOURCE_DIR}/opm/core/grid/cpgpreprocess/mxgrdecl.c
	${PROJECT_SOURCE_DIR}/opm/core/grid/cpgpreprocess/processgrid.c
	${PROJECT_SOURCE_DIR}/opm/core/utility/parameters/tinyxml/xmltest.cpp
	)

# remove inline TinyXML if a system version was found
if (TinyXML_FOUND)
	file (GLOB_RECURSE _inline_tinyxml "opm/core/utility/parameters/tinyxml/*")
	foreach (_file IN LISTS _inline_tinyxml)
		list (REMOVE_ITEM opmcore_SOURCES ${_file})
	endforeach (_file)
endif (TinyXML_FOUND)

# anyhow remove it from the header list (so it doesn't get installed)
list (REMOVE_ITEM opmcore_HEADERS "opm/core/utility/parameters/tinyxml/tinystr.h")
list (REMOVE_ITEM opmcore_HEADERS "opm/core/utility/parameters/tinyxml/tinyxml.h")

# HAVE_ERT is used as an #ifdef, not as an #if in the source code, if it
# is not true, then it should be unset altogether
if (NOT HAVE_ERT)
	set (HAVE_ERT)
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

# some CMake properties do not do list expansion
string (REPLACE ";" " " opm-core_LINKER_FLAGS_STR "${opm-core_LINKER_FLAGS}")

# default to building a shared library, but let user override
if (DEFINED BUILD_SHARED_LIBS)
	set (_shared_def ${BUILD_SHARED_LIBS})
else (DEFINED BUILD_SHARED_LIBS)
	set (_shared_def ON)
endif (DEFINED BUILD_SHARED_LIBS)
option (BUILD_opm-core_SHARED "Build opm-core as a shared library" ${_shared_def})
if (BUILD_opm-core_SHARED)
	set (opm-core_LIBRARY_TYPE SHARED)
else (BUILD_opm-core_SHARED)
	set (opm-core_LIBRARY_TYPE STATIC)
endif (BUILD_opm-core_SHARED)

# create this library
include_directories (${opm-core_INCLUDE_DIRS})
link_directories (${opm-core_LIBRARY_DIRS})
add_definitions (${opm-core_DEFINITIONS})
add_library (opmcore ${opm-core_LIBRARY_TYPE} ${opmcore_SOURCES})
set (opm-core_VERSION "${opm-core_VERSION_MAJOR}.${opm-core_VERSION_MINOR}")
set_target_properties (opmcore PROPERTIES
	SOVERSION ${opm-core_VERSION_MAJOR}
	VERSION ${opm-core_VERSION}
	LINK_FLAGS "${opm-core_LINKER_FLAGS_STR}"
	)
target_link_libraries (opmcore ${opm-core_LIBRARIES})

# queue this executable to be stripped
strip_debug_symbols (opmcore opm-core_DEBUG)

# pre-compile common headers; this is setup *after* the library to pick
# up extra options set there
option (PRECOMPILE_HEADERS "Precompile common headers for speed." ON)
mark_as_advanced (PRECOMPILE_HEADERS)
if (PRECOMPILE_HEADERS)
	get_target_property (_type opmcore TYPE)
	precompile_header (CXX ${_type}
		HEADER "${opmcore_PRECOMP_CXX_HEADER}"
		TARGET opmcore_CXX_pch
		FLAGS  opmcore_PRECOMP_CXX_FLAGS
		)
	# must set property on source files instead of entire target, because
	# it only applies to C++ modules (and cannot be used for C)
	set_source_files_properties (${opmcore_CXX_SOURCES} PROPERTIES
		OBJECT_DEPENDS "${opmcore_CXX_pch}"
		COMPILE_FLAGS  "${opmcore_PRECOMP_CXX_FLAGS}"
		)
	message (STATUS "Precompiled headers: ${opmcore_CXX_pch}")
else (PRECOMPILE_HEADERS)
	message (STATUS "Precompiled headers: disabled")
endif (PRECOMPILE_HEADERS)

### installation ###
foreach (_hdr IN LISTS opmcore_HEADERS)
	get_filename_component (_dir ${_hdr} PATH)
	file (RELATIVE_PATH _rel_dir "${PROJECT_SOURCE_DIR}" "${_dir}")
	install (
		FILES ${_hdr}
		DESTINATION include/${_rel_dir}
		)
endforeach (_hdr)
install (
	TARGETS opmcore
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	)
# only /usr/lib/debug seems to be searched for debug info; if we have
# write access to that directory (package installation), then default
# to use it; otherwise put the debug files together with the library
# (local installation). everything can be overridden by the option.
if (CMAKE_INSTALL_PREFIX STREQUAL "/usr")
	set (_sys_dbg_def ON)
else (CMAKE_INSTALL_PREFIX STREQUAL "/usr")
	set (_sys_dbg_def OFF)
endif (CMAKE_INSTALL_PREFIX STREQUAL "/usr")
option (SYSTEM_DEBUG "Put .debug files in GDB debug file directory" ${_sys_dbg_def})
set (DEBUG_FILE_DIRECTORY /usr/lib/debug CACHE LOCATION "GDB debug file directory")
mark_as_advanced (DEBUG_FILE_DIRECTORY)
if (SYSTEM_DEBUG)
	set (_dbg_prefix "${DEBUG_FILE_DIRECTORY}/")
else (SYSTEM_DEBUG)
	set (_dbg_prefix "")
endif (SYSTEM_DEBUG)
install (
	FILES ${PROJECT_BINARY_DIR}/${opm-core_DEBUG}
	DESTINATION ${_dbg_prefix}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
	)
install (
	FILES ${PROJECT_BINARY_DIR}/dune.module
	DESTINATION ${CMAKE_INSTALL_LIBDIR_NOARCH}/dunecontrol/opm-core
	)
message (STATUS "This build defaults to installing in ${CMAKE_INSTALL_PREFIX}")

# we need to know the name of the library which is generated
get_target_property (opm-core_LIBRARY opmcore LOCATION)

# installation of CMake modules to help user programs locate the library
include (OpmProject)
opm_cmake_config (opm-core)

### tutorials ###

# enumerate tutorials in project
file (GLOB tutorial_SOURCES "tutorials/tutorial[0-9].cpp")

# add each tutorial as a target of its own
add_custom_target (tutorials)
foreach (tutorial_FILE IN LISTS tutorial_SOURCES)
	# get the name of the executable based on the source code file
	get_filename_component (tutorial_NAME "${tutorial_FILE}" NAME_WE)

	# setup an executable, linking to the main library
	add_executable (${tutorial_NAME} EXCLUDE_FROM_ALL ${tutorial_FILE})
	set_target_properties (${tutorial_NAME} PROPERTIES
		LINK_FLAGS "${opm-core_LINKER_FLAGS_STR}"
		)
	target_link_libraries (${tutorial_NAME} opmcore ${opm-core_LIBRARIES})
	strip_debug_symbols (${tutorial_NAME})

	# add to list of tutorials to build
	add_dependencies (tutorials ${tutorial_NAME})
endforeach (tutorial_FILE)

### test programs ###
enable_testing ()
include (CTest)

# enumerate all testing programs
file (GLOB_RECURSE tests_SOURCES "tests/test_*.cpp")
file (GLOB_RECURSE not_tests_SOURCES "tests/not-unit/test_*.cpp")
list (REMOVE_ITEM tests_SOURCES ${not_tests_SOURCES})

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

# if ever huge test datafiles are necessary, then change this
# into "create_symlink" (on UNIX only, apparently)
set (make_avail "copy")

# provide datafiles as inputs for the tests, by copying them
# to a tests/ directory in the output tree (if different)
set (test_INPUT_FILES)
file (GLOB test_DATA "tests/*.xml")
foreach (input_datafile IN LISTS test_DATA)
	file (RELATIVE_PATH rel_datafile "${PROJECT_SOURCE_DIR}" ${input_datafile})
	set (output_datafile "${PROJECT_BINARY_DIR}/${rel_datafile}")
	if (NOT output_datafile STREQUAL input_datafile)
		add_custom_command (
			OUTPUT ${output_datafile}
			COMMAND ${CMAKE_COMMAND}
			ARGS -E ${make_avail} ${input_datafile} ${output_datafile}
			DEPENDS ${input_datafile}
			VERBATIM
			)
	endif (NOT output_datafile STREQUAL input_datafile)
	list (APPEND test_INPUT_FILES "${output_datafile}")
endforeach (input_datafile)
add_custom_target (datafiles
	DEPENDS ${test_INPUT_FILES}
	COMMENT "Making test data available in output tree"
	)

# compile each of these separately
add_custom_target (tests DEPENDS datafiles)
foreach (test_FILE IN LISTS tests_SOURCES)
	get_filename_component (test_NAME "${test_FILE}" NAME_WE)
	add_executable (${test_NAME} EXCLUDE_FROM_ALL ${test_FILE})
	add_dependencies (tests ${test_NAME})
	set_target_properties (${test_NAME} PROPERTIES
		LINK_FLAGS "${opm-core_LINKER_FLAGS_STR}"
		)
	target_link_libraries (${test_NAME} opmcore ${opm-core_LIBRARIES})
	strip_debug_symbols (${test_NAME})
	string (REGEX REPLACE "^test_([^/]*)$" "\\1" test_FANCY "${test_NAME}")
	get_target_property (test_LOC ${test_NAME} LOCATION)
	add_test (${test_FANCY} ${test_LOC})
	# run the test in the directory where the data files are
	set_tests_properties (${test_FANCY} PROPERTIES
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
		)
endforeach (test_FILE)

# ideally this should run the tests, but for now we only test if they
# compile without errors
add_custom_target (check
	COMMAND ${CMAKE_CTEST_COMMAND}
	DEPENDS tests
	COMMENT Checking if library is functional
	VERBATIM
	)

### documentation ###
configure_file (
	${PROJECT_SOURCE_DIR}/Doxyfile.in
	${PROJECT_BINARY_DIR}/Doxyfile
	@ONLY
	)
find_package (Doxygen)
if (DOXYGEN_FOUND)
	add_custom_target (doc
		COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
		SOURCES ${PROJECT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen"
		VERBATIM
		)
	set (_formats html)
	foreach (format IN LISTS _formats)
		string (TOUPPER ${format} FORMAT)
		install (
			DIRECTORY ${PROJECT_BINARY_DIR}/Documentation/${format}
			DESTINATION share/doc/opm-core/
			COMPONENT ${format}
			OPTIONAL
			)
		# target to install just HTML documentation
		add_custom_target (install-${format}
			COMMAND ${CMAKE_COMMAND} -DCOMPONENT=${format} -P cmake_install.cmake
			COMMENT Installing ${FORMAT} documentation
			VERBATIM
			)
		# since the documentation is optional, it is not automatically built
		add_dependencies (install-${format} doc)
	endforeach (format)
endif (DOXYGEN_FOUND)

### clean in-source builds ###
set (DISTCLEAN_FILES
	CMakeCache.txt
	cmake_install.cmake
	Makefile
  config.h
  opm-core-config.cmake
	opm-core-config-version.cmake
  opm-core-install.cmake
  Doxyfile
	)
set (DISTCLEAN_DIRS
	CMakeFiles/
	)
add_custom_target (distclean
	COMMAND ${CMAKE_COMMAND} -E remove ${DISTCLEAN_FILES}
	COMMAND ${CMAKE_COMMAND} -E remove_directory ${DISTCLEAN_DIRS}
# cannot depend on clean because it is only defined in the master Makefile
# not in CMakeFiles/Makefile where this target will end up
#	DEPENDS clean
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	COMMENT Removing CMake-generated files
	VERBATIM
	)

# smart wrapper that auto-parallelizes builds
file (COPY
	GNUmakefile
	DESTINATION ${PROJECT_BINARY_DIR}
	)
