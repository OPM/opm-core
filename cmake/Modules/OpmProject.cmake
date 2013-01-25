# - Helper routines for opm-core like projects

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
  configure_file (
	${PROJECT_SOURCE_DIR}/dune.module.in
	${PROJECT_BINARY_DIR}/dune.module
	@ONLY
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
endfunction (opm_cmake_config name)
