# - Provide C wrappers for Fortran code
#
# Synopsis:
#	define_fc_func (APPEND config.h)

function (define_fc_func verb file)
  # check that we are being called correctly
  if (NOT (("${verb}" STREQUAL "APPEND") OR
		   ("${verb}" STREQUAL "WRITE")))
	message (FATAL_ERROR
	  "Unknown verb \"${verb}\" passed as first argument."
	  )
  endif (NOT (("${verb}" STREQUAL "APPEND") OR
		      ("${verb}" STREQUAL "WRITE")))

  # only do something if we actually have some components which requires
  # the interaction -- don't load the Fortran compiler just to write
  # this macro (which apparently nobody uses then)
  if (CMAKE_C_COMPILER_LOADED AND CMAKE_Fortran_COMPILER_LOADED)
	# get a temporary file
	set (_tmp_hdr ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/config_f.h)

	# write a small config file that contains the proper convention for
	# calling Fortran from C
	include (FortranCInterface)
	fortrancinterface_header (${_tmp_hdr})

	# read the definition back in from the file
	file (STRINGS
	  ${_tmp_hdr}
	  _str
	  REGEX "^#define FortranCInterface_GLOBAL\\(name,NAME\\) .*$"
	  )

	# massage it to look like the one AC_FC_WRAPPERS provide
	string (REPLACE "FortranCInterface_GLOBAL" "FC_FUNC" _str ${_str})

	# write this definition to the end of our own configuration file
	file (${verb} ${file}
	  "\n/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */\n"
      ${_str}
	  "\n"
	  )
  else (CMAKE_C_COMPILER_LOADED AND CMAKE_Fortran_COMPILER_LOADED)
	message (STATUS "Fortran/C interface not activated")
  endif (CMAKE_C_COMPILER_LOADED AND CMAKE_Fortran_COMPILER_LOADED)
endfunction (define_fc_func)
