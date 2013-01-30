# - Helper routines for opm-core like projects

# convert a list back to a command-line string
function (unseparate_args var_name prefix value)
  separate_arguments (value)
  foreach (item IN LISTS value)
	set (prefixed_item "${prefix}${item}")
	if (${var_name})
	  set (${var_name} "${${var_name}} ${prefixed_item}")
	else (${var_name})
	  set (${var_name} "${prefixed_item}")
	endif (${var_name})
  endforeach (item)
  set (${var_name} "${${var_name}}" PARENT_SCOPE)
endfunction (unseparate_args var_name prefix value)

# wrapper to set variables in pkg-config file
function (configure_pc_file name source dest libdir includedir)
  # escape set of standard strings
  unseparate_args (includes "-I" "${${name}_INCLUDE_DIRS}")
  unseparate_args (libs "-l" "${${name}_LIBRARIES}")
  unseparate_args (defs "" "${${name}_DEFINITIONS}")

  # necessary to make these variables visible to configure_file
  set (target "${${name}_LIBRARY}")
  set (major "${${name}_VERSION_MAJOR}")
  set (minor "${${name}_VERSION_MINOR}")

  configure_file (${source} ${dest} @ONLY)
endfunction (configure_pc_file name source dist libdir includedir)

# installation of CMake modules to help user programs locate the library
function (opm_cmake_config name)
  # write configuration file to locate library
  configure_file (
	${PROJECT_SOURCE_DIR}/${name}-config.cmake.in
	${PROJECT_BINARY_DIR}/${name}-config.cmake
	@ONLY
	)
  configure_file (
	${PROJECT_SOURCE_DIR}/${name}-config-version.cmake.in
	${PROJECT_BINARY_DIR}/${name}-config-version.cmake
	@ONLY
	)
  configure_vars (
	FILE CMAKE "${PROJECT_BINARY_DIR}/${name}-config.cmake"
	APPEND "${${name}_CONFIG_VARS}"
	)

  # config-mode .pc file; use this to find the build tree
  configure_pc_file (
	${name}
	${PROJECT_SOURCE_DIR}/${name}.pc.in
	${PROJECT_BINARY_DIR}/${name}.pc
	${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
	${PROJECT_SOURCE_DIR}
	)

  # replace the build directory with the target directory in the
  # variables that contains build paths
  string (REPLACE
	"${PROJECT_SOURCE_DIR}"
	"${CMAKE_INSTALL_PREFIX}/include"
	${name}_INCLUDE_DIRS
	"${${name}_INCLUDE_DIRS}"
	)
  string (REPLACE
	"${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
	"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
	${name}_LIBRARY
	"${${name}_LIBRARY}"
	)
  set (CMAKE_LIBRARY_OUTPUT_DIRECTORY
	"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
	)
  # create a config mode file which targets the install directory instead
  # of the build directory (using the same input template)
  configure_file (
	${PROJECT_SOURCE_DIR}/${name}-config.cmake.in
	${PROJECT_BINARY_DIR}/${name}-install.cmake
	@ONLY
	)
  configure_vars (
	FILE CMAKE "${PROJECT_BINARY_DIR}/${name}-install.cmake"
	APPEND "${${name}_CONFIG_VARS}"
	)
  # this file gets copied to the final installation directory
  install (
	FILES ${PROJECT_BINARY_DIR}/${name}-install.cmake
	DESTINATION share/cmake/${name}
	RENAME ${name}-config.cmake
	)
  # assume that there exists a version file already
  install (
	FILES ${PROJECT_BINARY_DIR}/${name}-config-version.cmake
	DESTINATION share/cmake/${name}
	)

  # find-mode .pc file; use this to locate system installation
  configure_pc_file (
	${name}
	${PROJECT_SOURCE_DIR}/${name}.pc.in
	${PROJECT_BINARY_DIR}/${name}-install.pc
	${CMAKE_INSTALL_LIBDIR}
	${CMAKE_INSTALL_PREFIX}/include
	)

  # put this in the right system location; assume that we have binaries
  install (
	FILES ${PROJECT_BINARY_DIR}/${name}-install.pc
	DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/pkgconfig/
	RENAME ${name}.pc
	)
endfunction (opm_cmake_config name)
