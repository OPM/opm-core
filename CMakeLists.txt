# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:
cmake_minimum_required (VERSION 2.8)

# key information about the library
set (project "opm-core")
set (${project}_NAME "${project}")
set (${project}_DESCRIPTION "Open Porous Media Initiative Core Library")
set (${project}_DIR "opm")
set (${project}_VERSION_MAJOR 1)
set (${project}_VERSION_MINOR 0)
set (doxy_dir "Documentation")

# additional search modules
set (${project}_MODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
list (APPEND CMAKE_MODULE_PATH ${${project}_MODULE_DIR})

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include (${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
	# Algebraic Multigrid must be compiled together with our program;
	# if it is not available, then remove our corresponding component
	find_package (AGMG)
	if (AGMG_FOUND)
		list (APPEND ${project}_SOURCES ${AGMG_SOURCES})
	endif (AGMG_FOUND)

	# these solvers are only compiled in if their dependency is found
	if (NOT AGMG_FOUND)
		list (REMOVE_ITEM opm-core_SOURCES
			${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/LinearSolverAGMG.cpp
			)
	endif (NOT AGMG_FOUND)
	if (NOT dune-istl_FOUND)
		list (REMOVE_ITEM opm-core_SOURCES
			${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/LinearSolverIstl.cpp
			)
	endif (NOT dune-istl_FOUND)
	if (NOT SuiteSparse_FOUND)
		list (REMOVE_ITEM opm-core_SOURCES
			${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/call_umfpack.c
			${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/LinearSolverUmfpack.cpp
			)
		list (REMOVE_ITEM examples_SOURCES
			${PROJECT_SOURCE_DIR}/tutorials/tutorial2.cpp
			${PROJECT_SOURCE_DIR}/tutorials/tutorial3.cpp
			${PROJECT_SOURCE_DIR}/tutorials/tutorial4.cpp
			)
	endif (NOT SuiteSparse_FOUND)

	# we are not supposed to include the TinyXML test prog. regardless
	list (REMOVE_ITEM opm-core_SOURCES
		${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/utility/parameters/tinyxml/xmltest.cpp
		)

	# remove inline TinyXML if a system version was found
	if (TinyXML_FOUND)
		file (GLOB_RECURSE _inline_tinyxml "${opm-core_DIR}/core/utility/parameters/tinyxml/*")
		foreach (_file IN LISTS _inline_tinyxml)
			list (REMOVE_ITEM opm-core_SOURCES ${_file})
		endforeach (_file)
	endif (TinyXML_FOUND)

	# anyhow remove it from the header list (so it doesn't get installed)
	list (REMOVE_ITEM opm-core_HEADERS "${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/utility/parameters/tinyxml/tinystr.h")
	list (REMOVE_ITEM opm-core_HEADERS "${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/utility/parameters/tinyxml/tinyxml.h")

	# HAVE_ERT is used as an #ifdef, not as an #if in the source code, if it
	# is not true, then it should be unset altogether
	if (NOT HAVE_ERT)
		set (HAVE_ERT)
		list (REMOVE_ITEM examples_SOURCES
			${PROJECT_SOURCE_DIR}/examples/import_rewrite.cpp
			)
	endif (NOT HAVE_ERT)
endmacro (sources_hook)

macro (fortran_hook)
	# only include Fortran support if AGMG sources are available
	set (${project}_FORTRAN_IF HAVE_AGMG)
endmacro (fortran_hook)

macro (tests_hook)
	cond_disable_test ("AGMG")
	cond_disable_test ("ERT")
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)
