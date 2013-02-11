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
  option (BUILD_${${opm}_NAME}_SHARED "Build ${${opm}_NAME} as a shared library" ${_shared_def})
  if (BUILD_${${opm}_NAME}_SHARED)
	set (${opm}_LIBRARY_TYPE SHARED)
  else (BUILD_${${opm}_NAME}_SHARED)
	set (${opm}_LIBRARY_TYPE STATIC)
  endif (BUILD_${${opm}_NAME}_SHARED)

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
endmacro (opm_defaults opm)
