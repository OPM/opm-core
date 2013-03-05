# - Default settings for the build

include (UseCompVer)

macro (opm_defaults opm)
  # build debug by default
  if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Debug")
  endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

  # default to building a static library, but let user override
  if (DEFINED BUILD_SHARED_LIBS)
	set (_shared_def ${BUILD_SHARED_LIBS})
  else (DEFINED BUILD_SHARED_LIBS)
	set (_shared_def OFF)
  endif (DEFINED BUILD_SHARED_LIBS)
  string (TOUPPER "${${opm}_NAME}" opm_UPPER)
  string (REPLACE "-" "_" opm_UPPER "${opm_UPPER}")
  option (BUILD_${opm_UPPER}_SHARED "Build ${${opm}_NAME} as a shared library" ${_shared_def})
  if (BUILD_${opm_UPPER}_SHARED)
	set (${opm}_LIBRARY_TYPE SHARED)
  else (BUILD_${opm_UPPER}_SHARED)
	set (${opm}_LIBRARY_TYPE STATIC)
  endif (BUILD_${opm_UPPER}_SHARED)

  # precompile standard headers to speed up compilation
  # unfortunately, this functionality is buggy and tends to segfault at
  # least up to version 4.7.2, so it should be disabled by default there
  get_gcc_version (CXX GCC_VERSION)
  if (GCC_VERSION VERSION_LESS "4.7.2")
	set (_precomp_def OFF)
  else (GCC_VERSION VERSION_LESS "4.7.2")
	set (_precomp_def ON)
  endif (GCC_VERSION VERSION_LESS "4.7.2")
  option (PRECOMPILE_HEADERS "Precompile common headers for speed." ${_precomp_def})
  mark_as_advanced (PRECOMPILE_HEADERS)
  if (NOT PRECOMPILE_HEADERS)
	message (STATUS "Precompiled headers: disabled")
  endif(NOT PRECOMPILE_HEADERS)

  # if we are on a system where CMake 2.6 is the default (Hi RHEL 6!),
  # the configuration files for Boost will trip up the library paths
  # (look for /usr/lib64/lib64/ in the log) when used with FindBoost
  # module bundled with CMake 2.8. this can be circumvented by turning
  # off config mode probing if we have not explicitly specified a
  # directory to look for it. for more details, see
  # <http://stackoverflow.com/questions/9948375/cmake-find-package-succeeds-but-returns-wrong-path>
  if (NOT BOOST_ROOT)
	set (Boost_NO_BOOST_CMAKE ON)
  endif (NOT BOOST_ROOT)
endmacro (opm_defaults opm)
