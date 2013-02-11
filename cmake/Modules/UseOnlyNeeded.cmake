# - Use only needed imports from libraries
#
# Add the -Wl,--as-needed flag to the default linker flags on Linux
# in order to get only the minimal set of dependencies.

function (prepend var_name value)
  if (${var_name})
	set (${var_name} "${value} ${${var_name}}" PARENT_SCOPE)
  else (${var_name})
	set (${var_name} "${value}")
  endif (${var_name})
endfunction (prepend var_name value)

if (CMAKE_CXX_PLATFORM_ID STREQUAL "Linux")
  prepend (CMAKE_EXE_LINKER_FLAGS "-Wl,--as-needed")
  prepend (CMAKE_MODULE_LINKER_FLAGS "-Wl,--as-needed")
  prepend (CMAKE_SHARED_LINKER_FLAGS "-Wl,--as-needed")
endif (CMAKE_CXX_PLATFORM_ID STREQUAL "Linux")
