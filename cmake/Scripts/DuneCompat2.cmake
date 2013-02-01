# - Emulate a rule to patch the Makefile, adding a line to the source
# tree and write a marker file indicating it is done.

set (base_dir ".")
set (marker_file "${base_dir}/CMakeFiles/marker")
set (makefile "${base_dir}/Makefile")

# if the Makefile has changed, then update it and touch the marker so that
# the marker is newer and won't update the Makefile again
if ("${makefile}" IS_NEWER_THAN "${marker_file}")
  file (APPEND "${makefile}"
	"abs_top_srcdir = ${CMAKE_HOME_DIRECTORY}\n"
	)
  execute_process (COMMAND
	${CMAKE_COMMAND} -E touch "${marker_file}"
	)
endif ("${makefile}" IS_NEWER_THAN "${marker_file}")
	
