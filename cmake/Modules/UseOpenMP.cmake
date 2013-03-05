# - Use OpenMP features
#
# Synopsis:
#
#	find_openmp (module)
#
# where:
#
#	module            Name of the module to which OpenMP libraries
#	                  etc. should be added, e.g. "opm-core".
#
# Note: Compiler flags are always added globally, to avoid ABI
# incompatibility problems.
#
# It is assumed that the following variables are available
#
#	${module}_QUIET      Verbosity level of the parent's find module
#	${module}_LIBRARIES  List of libraries to which OpenMP will be added
#
# Example:
#	find_openmp (opm-core)
#	remove_dup_deps (opm-core)

include (AddOptions)
macro (find_openmp opm)
  # enabling OpenMP is supposedly enough to make the compiler link with
  # the appropriate libraries
  find_package (OpenMP ${${opm}_QUIET})
  list (APPEND ${opm}_LIBRARIES ${OpenMP_LIBRARIES})
  if (OPENMP_FOUND)
	add_options (C ALL_BUILDS "${OpenMP_C_FLAGS}")
	add_options (CXX ALL_BUILDS "${OpenMP_CXX_FLAGS}")
  endif (OPENMP_FOUND)

  # threading library (search for this *after* OpenMP
  set (CMAKE_THREAD_PREFER_PTHREAD TRUE)
  find_package (Threads ${${opm}_QUIET})
  if (CMAKE_USE_PTHREADS_INIT)
	list (APPEND ${opm}_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
  endif (CMAKE_USE_PTHREADS_INIT)
endmacro (find_openmp opm)
