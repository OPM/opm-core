# - Build satellites that are dependent of main library

# Synopsis:
#	opm_data (satellite target files)
#
# provides these output variables:
#
#	${satellite_INPUT_FILES}   List of all files that are copied
#	${satellite}_DATAFILES     Name of target which copies these files
#
# Example:
#
#	opm_data (test datafiles "tests/*.xml")
#
macro (opm_data satellite target files)
  # if ever huge test datafiles are necessary, then change this
  # into "create_symlink" (on UNIX only, apparently)
  set (make_avail "copy")

  # provide datafiles as inputs for the tests, by copying them
  # to a tests/ directory in the output tree (if different)
  set (${satellite}_INPUT_FILES)
  file (GLOB ${satellite}_DATA "tests/*.xml")
  if (NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
	foreach (input_datafile IN LISTS ${satellite}_DATA)
	  file (RELATIVE_PATH rel_datafile "${PROJECT_SOURCE_DIR}" ${input_datafile})
	  set (output_datafile "${PROJECT_BINARY_DIR}/${rel_datafile}")
	  add_custom_command (
		OUTPUT ${output_datafile}
		COMMAND ${CMAKE_COMMAND}
		ARGS -E ${make_avail} ${input_datafile} ${output_datafile}
		DEPENDS ${input_datafile}
		VERBATIM
		)
	  list (APPEND ${satellite}_INPUT_FILES "${output_datafile}")
	endforeach (input_datafile)
  endif(NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)

  # setup a target which does all the copying
  set (${satellite}_DATAFILES "${target}")
  add_custom_target (${${satellite}_DATAFILES}
	DEPENDS ${${satellite}_INPUT_FILES}
	COMMENT "Making test data available in output tree"
	)
endmacro (opm_data satellite target files)
