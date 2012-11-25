# - Create config.h based on a list of variables
#
# Synopsis:
#   configure_vars (FILE filename verb varlist)
# where
#   filename      Full path (including name) of config.h
#   verb          WRITE or APPEND if truncating or not
#   varlist       List of variable names that has been defined
#
# In addition, this function will define HAVE_CONFIG_H for the
# following compilations, (only) if the filename is "config.h".
#
# Example:
#   list (APPEND FOO_CONFIG_VARS
#     "/* bar library */"
#     "HAVE_BAR"
#     "HAVE_BAR_VERSION_2"
#     )
#   configure_vars (
#     FILE  ${PROJECT_BINARY_DIR}/config.h
#     WRITE ${FOO_CONFIG_VARS}
#     )

# Copyright (C) 2012 Uni Research AS
# This file is licensed under the GNU General Public License v3.0

function (configure_vars obj filename verb)
  # this is just to make the syntax look like the build-in commands
  message (STATUS "Writing config file \"${filename}\"...")
  if (NOT ("${obj}" STREQUAL "FILE" AND
		(("${verb}" STREQUAL "WRITE") OR ("${verb}" STREQUAL "APPEND"))))
	message (FATAL_ERROR "Syntax error in argument list")
  endif (NOT ("${obj}" STREQUAL "FILE" AND
	  (("${verb}" STREQUAL "WRITE") OR ("${verb}" STREQUAL "APPEND"))))
  
  # truncate the file if the verb was "WRITE"
  if (verb STREQUAL "WRITE")
	file (WRITE "${filename}" "")
  endif (verb STREQUAL "WRITE")
  
  # whenever we use this, we also signal to the header files that we
  # have "config.h". add this before any other files (known till now)
  # to avoid confusion from other configuration files.
  get_filename_component (_config_path "${filename}" PATH)
  get_filename_component (_config_file "${filename}" NAME)
  if ("${_config_file}" STREQUAL "config.h")
	add_definitions (-DHAVE_CONFIG_H=1)
	include_directories (BEFORE "${_config_path}")
  endif ("${_config_file}" STREQUAL "config.h")
  
  # only write the current value of each variable once
  list (APPEND _args ${ARGN})
  list (REMOVE_DUPLICATES _args)
  
  # process each variable
  set (_prev_verbatim TRUE)
  foreach (_var IN LISTS _args)

	# massage the name to remove source code formatting
	string (REGEX REPLACE "^[\\n\\t\\ ]+" "" _var "${_var}")
	string (REGEX REPLACE "[\\n\\t\\ ]+$" "" _var "${_var}")
	
	# if the name of a variable has the syntax of a comments, write it
	# verbatim to the file; this can be used to create headings
	if ("${_var}" MATCHES "^/[/*]")
	  if (NOT _prev_verbatim)
		file (APPEND "${filename}" "\n")
	  endif (NOT _prev_verbatim)
	  file (APPEND "${filename}" "${_var}\n")
	  set (_prev_verbatim TRUE)
	  
	else ("${_var}" MATCHES "^/[/*]")
	  
	  # check for empty variable; variables that are explicitly set to false
	  # is not included in this clause
	  if ((NOT DEFINED ${_var}) OR ("${${_var}}" STREQUAL ""))
		file (APPEND "${filename}" "/* #undef ${_var} */\n")
	  else ((NOT DEFINED ${_var}) OR ("${${_var}}" STREQUAL ""))
		
		# integer variables (specifically 0 and 1) are written as they are,
		# whereas everything else (including version numbers, which could
		# be interpreted as floats) are quoted as strings
		if (${_var} MATCHES "[0-9]+")
		  file (APPEND "${filename}" "#define ${_var} ${${_var}}\n")
		else (${_var} MATCHES "[0-9]+")
		  file (APPEND "${filename}" "#define ${_var} \"${${_var}}\"\n")
		endif (${_var} MATCHES "[0-9]+")
		
	  endif ((NOT DEFINED ${_var}) OR ("${${_var}}" STREQUAL ""))
	  set (_prev_verbatim FALSE)
	endif ("${_var}" MATCHES "^/[/*]")
  endforeach(_var)
endfunction (configure_vars)
