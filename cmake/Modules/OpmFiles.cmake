# - Identify source code

macro (opm_out_dirs)
  # put libraries in lib/ (no multi-arch support in build tree)
  set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
  set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
  set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
  set (CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/CMakeFiles")
endmacro (opm_out_dirs)

macro (opm_sources opm)
  # find all the source code (note that these variables have name after
  # the target library and not the project). the documentation recommends
  # against using globs to enumerate source code, but if you list the
  # files explicitly you'll change the build files every time you add to
  # the project as well as having to rebuild completely anyway.
  file (GLOB_RECURSE ${opm}_C_SOURCES "${${opm}_DIR}/[^.]*.c")
  file (GLOB_RECURSE ${opm}_CXX_SOURCES "${${opm}_DIR}/[^.]*.cpp")
  file (GLOB_RECURSE ${opm}_C_HEADERS "${${opm}_DIR}/[^.]*.h")
  file (GLOB_RECURSE ${opm}_CXX_HEADERS "${${opm}_DIR}/[^.]*.hpp")

  # remove pre-compile headers from output list
  set (${opm}_PRECOMP_CXX_HEADER "${${opm}_DIR}/${${opm}_NAME}-pch.hpp")
  list (REMOVE_ITEM ${opm}_CXX_HEADERS
	${PROJECT_SOURCE_DIR}/${${opm}_PRECOMP_CXX_HEADER}
	)

  # merge both languages into one compilation/installation
  set (${opm}_SOURCES ${${opm}_C_SOURCES} ${${opm}_CXX_SOURCES})
  set (${opm}_HEADERS ${${opm}_C_HEADERS} ${${opm}_CXX_HEADERS})
endmacro (opm_sources opm)

macro (opm_find_tests)
  # every C++ program prefixed with test_ under tests/ directory should
  # be automatically set up as a test
  set (tests_DIR "tests")
  file (GLOB_RECURSE tests_SOURCES "${tests_DIR}/test_*.cpp")
  file (GLOB_RECURSE not_tests_SOURCES "${tests_DIR}/not-unit/test_*.cpp")
  list (REMOVE_ITEM tests_SOURCES ${not_tests_SOURCES})
endmacro (opm_find_tests)

macro (opm_find_tutorials)
  # enumerate tutorials in project
  set (tutorial_DIR "tutorials")
  file (GLOB tutorial_SOURCES "${tutorial_DIR}/tutorial[0-9].cpp")
endmacro (opm_find_tutorials)
