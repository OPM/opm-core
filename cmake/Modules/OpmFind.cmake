# - Generic inclusion of packages
#
# Synopsis:
#
#	find_and_append_package (name args)
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
# - Generic inclusion of a list of packages
#
# Synopsis:
#
#	find_and_append_package_list (args)
#
# where
#
#	args          List of package strings. Each string must be quoted if
#	              it contains more than one word.
#
# Example:
#
#	find_and_append_package_list (
#		"Boost COMPONENTS filesystem REQUIRED"
#		SUPERLU
#	)
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

# find a list of dependencies, adding each one of them
macro (find_and_append_package_list_to prefix)
  # setting and separating is necessary to work around apparent bugs
  # in CMake's parser (sic)
  set (_deps ${ARGN})
  foreach (_dep IN LISTS _deps)
	separate_arguments (_args UNIX_COMMAND ${_dep})
	find_and_append_package_to (${prefix} ${_args})
  endforeach (_dep)
endmacro (find_and_append_package_list_to prefix)

# convenience method to supply the project name as prefix
macro (find_and_append_package_list)
  find_and_append_package_list_to (${CMAKE_PROJECT_NAME} ${ARGN})
endmacro (find_and_append_package_list)

# libraries should always be trimmed from the beginning, so that also
# missing functions in those later in the list will be resolved
macro (remove_duplicate_libraries module)
  if (DEFINED ${module}_LIBRARIES)
	list (REVERSE ${module}_LIBRARIES)
	list (REMOVE_DUPLICATES ${module}_LIBRARIES)
	list (REVERSE ${module}_LIBRARIES)
  endif (DEFINED ${module}_LIBRARIES)
endmacro (remove_duplicate_libraries module)
