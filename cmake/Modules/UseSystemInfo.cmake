# - Print CMake and OS distribution version
#
function (system_info)
  message (STATUS "CMake version: ${CMAKE_VERSION}")
  if (CMAKE_SYSTEM MATCHES "Linux")
	read_lsb (ID INTO DISTRIB_ID)
	read_lsb (RELEASE INTO DISTRIB_RELEASE)
	read_lsb (CODENAME INTO DISTRIB_CODENAME)
	message (STATUS "Linux distribution: ${DISTRIB_ID} \"${DISTRIB_CODENAME}\" ${DISTRIB_RELEASE}")
  else (CMAKE_SYSTEM MATCHES "Linux")
	message (STATUS "Operating system: ${CMAKE_SYSTEM}")
  endif (CMAKE_SYSTEM MATCHES "Linux")
endfunction (system_info)

# read property from LSB information
function (read_lsb suffix INTO varname)
  file (STRINGS /etc/lsb-release _distrib
	REGEX "^DISTRIB_${suffix}="
	)
  string (REGEX REPLACE
	"^DISTRIB_${suffix}=\(.*\)" "\\1" ${varname} ${_distrib})
  set (${varname} "${${varname}}" PARENT_SCOPE)
endfunction (read_lsb suffix INTO varname)
