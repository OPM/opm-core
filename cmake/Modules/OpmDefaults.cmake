# - Default settings for the build

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
endmacro (opm_defaults opm)
