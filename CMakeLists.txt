# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:
cmake_minimum_required (VERSION 2.8)
project (opm-core)
enable_language (CXX)

# all public header files are together with the source
set (opm-core_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}")

# put libraries in lib/
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

# find all the source code
file (GLOB_RECURSE opm-core_SOURCES "opm/*.c" "opm/*.cpp")

# create this library
include_directories (${opm-core_INCLUDE_DIRS})
add_library (opmcore SHARED ${opm-core_SOURCES})
