# - Get compiler version

# probe the GCC version, returns empty string if GCC is not compiler
function (get_gcc_version language ver_name)
  if(CMAKE_${language}_COMPILER_ID STREQUAL GNU)  
	# exec_program is deprecated, but execute_process does't work :-(
	exec_program (${CMAKE_${language}_COMPILER}
	  ARGS ${CMAKE_${language}_COMPILER_ARG1} -dumpversion
	  OUTPUT_VARIABLE _version
	  )
	set (${ver_name} ${_version} PARENT_SCOPE)
  else (CMAKE_${language}_COMPILER_ID STREQUAL GNU)
	set (${ver_name} "" PARENT_SCOPE)
  endif (CMAKE_${language}_COMPILER_ID STREQUAL GNU)
endfunction (get_gcc_version ver_name)
