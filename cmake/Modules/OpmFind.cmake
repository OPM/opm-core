# - Generic inclusion of packages
#
# Synopsis:
#
#	find_and_append (name args)
#
# where
#
#	name          Name of the package, e.g. Boost
#   args          Other arguments, e.g. COMPONENTS, REQUIRED, QUIET etc.
#
# This macro will append the list of standard variables found by the
# package to this project's standard variables
#
########################################################################
#
# - Remove duplicate library declarations
#
# Synopsis:
#
#	remove_duplicate_libraries (module)
#
# where
#	module         Name of the module whose libraries should be pruned

# list of suffixes for all the project variables
set (_opm_proj_vars
  LINKER_FLAGS
  LIBRARIES
  DEFINITIONS
  INCLUDE_DIRS
  LIBRARY_DIRS
  CONFIG_VARS
  )

# ensure that they are at least the empty list after we're done
foreach (name IN LISTS _opm_proj_vars)
  if (NOT DEFINED ${CMAKE_PROJECT_NAME}_${name})
	set (${CMAKE_PROJECT_NAME}_${name} "")
  endif (NOT DEFINED ${CMAKE_PROJECT_NAME}_${name})
endforeach (name)

# insert this boilerplate whenever we are going to find a new package
macro (find_and_append_package_to prefix name)
  find_package (${name} ${ARGN})
  if (${name}_FOUND)
	foreach (var IN LISTS _opm_proj_vars)
	  if (DEFINED ${name}_${var})
		list (APPEND ${prefix}_${var} ${${name}_${var}})
	  endif (DEFINED ${name}_${var})
	endforeach (var)
  endif (${name}_FOUND)
endmacro (find_and_append_package_to prefix name)

# append to the list of variables associated with the project
macro (find_and_append_package name)
  find_and_append_package_to (${CMAKE_PROJECT_NAME} ${name} ${ARGN})
endmacro (find_and_append_package name)

# libraries should always be trimmed from the beginning, so that also
# missing functions in those later in the list will be resolved
macro (remove_duplicate_libraries module)
  if (DEFINED ${module}_LIBRARIES)
	list (REVERSE ${module}_LIBRARIES)
	list (REMOVE_DUPLICATES ${module}_LIBRARIES)
	list (REVERSE ${module}_LIBRARIES)
  endif (DEFINED ${module}_LIBRARIES)
endmacro (remove_duplicate_libraries module)
