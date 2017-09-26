# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################
 
# Mandatory call to project
project(opm-core C CXX)

cmake_minimum_required (VERSION 2.8)

# additional search modules
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common}")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

include(OpmInit)

# not the same location as most of the other projects; this hook overrides
macro (dir_hook)
	set (doxy_dir "Documentation")
endmacro (dir_hook)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include ( ${project}-prereqs )

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
	opm_need_version_of ("dune-common")
  opm_need_version_of ("dune-istl")
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
	# these solvers are only compiled in if their dependency is found
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

	if (NOT PETSC_FOUND)
		list (REMOVE_ITEM opm-core_SOURCES
			${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/call_petsc.c
			${PROJECT_SOURCE_DIR}/${opm-core_DIR}/core/linalg/LinearSolverPetsc.cpp
			)
	endif (NOT PETSC_FOUND)
	if ((NOT MPI_FOUND) OR (NOT DUNE_ISTL_FOUND))
		list (REMOVE_ITEM tests_SOURCES
			${PROJECT_SOURCE_DIR}/tests/test_parallel_linearsolver.cpp
			${PROJECT_SOURCE_DIR}/tests/test_parallelistlinformation.cpp
			)
	endif ((NOT MPI_FOUND) OR (NOT DUNE_ISTL_FOUND))

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
endmacro (fortran_hook)

macro (files_hook)
endmacro (files_hook)

macro (tests_hook)
	cond_disable_test ("ERT")
endmacro (tests_hook)

macro (install_hook)
#	# copy the cmake macros to system when installing
#	install (
#		DIRECTORY cmake
#		DESTINATION share/opm
#		)
endmacro (install_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)
