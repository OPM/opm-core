# - Find the Ensemble-based Reservoir Tool (ERT)
#
# Set the cache variable ERT_ROOT to the install location of the ERT
# libraries and header files.
#
# If found, it sets these variables:
#
#	ERT_INCLUDE_DIRS      Header file directories
#	ERT_LIBRARIES         Archives and shared objects
#	ERT_CONFIG_VARS       Definitions that goes in config.h
#	ERT_LINKER_FLAGS      Options that must be passed to linker
#
# It will also add to CMAKE_C_FLAGS and CMAKE_CXX_FLAGS if necessary to
# link with the ERT libraries.

# variables to pass on to other packages
if (FIND_QUIETLY)
  set (ERT_QUIET "QUIET")
else (FIND_QUIETLY)
  set (ERT_QUIET "")
endif (FIND_QUIETLY)

# ERT doesn't have any config-mode file, so we need to specify the root
# directory in its own variable
find_path (ERT_INCLUDE_DIR
  NAMES "ecl_util.h"
  HINTS "${PROJECT_BINARY_DIR}/../ert"
  PATHS "${ERT_ROOT}"
  PATH_SUFFIXES "devel/libecl/src" "include"
  DOC "Path to ERT Eclipse library header files"
  )

# need all of these libraries
find_library (ERT_LIBRARY_ECL
  NAMES "ecl"
  PATHS "${ERT_ROOT}"
  PATH_SUFFIXES "devel/lib" "lib" "lib64" "lib32"
  DOC "Path to ERT Eclipse library archive/shared object files"
  )
find_library (ERT_LIBRARY_GEOMETRY
  NAMES "geometry"
  PATHS "${ERT_ROOT}"
  PATH_SUFFIXES "devel/lib" "lib" "lib64" "lib32"
  DOC "Path to ERT Geometry library archive/shared object files"
  )
find_library (ERT_LIBRARY_UTIL
  NAMES "ert_util"
  PATHS "${ERT_ROOT}"
  PATH_SUFFIXES "devel/lib" "lib" "lib64" "lib32"
  DOC "Path to ERT Utilities library archive/shared object files"
  )
# the "library" found here is actually a list of several files
list (APPEND ERT_LIBRARY
  ${ERT_LIBRARY_ECL}
  ${ERT_LIBRARY_GEOMETRY}
  ${ERT_LIBRARY_UTIL}
  )
list (APPEND ERT_LIBRARIES ${ERT_LIBRARY})
list (APPEND ERT_INCLUDE_DIRS ${ERT_INCLUDE_DIR})

# these system variables are probed for, and used in HEADER files (sic)
list (APPEND ERT_CONFIG_VARS
  HAVE_ISFINITE
  HAVE_GLOB
  HAVE_FORK
  HAVE_GETUID
  HAVE_LOCKF
  HAVE_OPENDIR
  HAVE_PROC
  HAVE_READLINKAT
  HAVE_SYMLINK
  HAVE_VA_COPY
  )
include (CheckSymbolExists)
include (CheckFunctionExists)
check_symbol_exists (isfinite math.h HAVE_ISFINITE)
check_function_exists (glob HAVE_GLOB)
check_function_exists (fork HAVE_FORK)
check_function_exists (getuid HAVE_GETUID)
check_function_exists (lockf HAVE_LOCKF)
check_function_exists (opendir HAVE_OPENDIR)
check_function_exists (readlinkat HAVE_READLINKAT)
check_function_exists (symlink HAVE_SYMLINK)
check_symbol_exists (va_copy stdarg.h HAVE_VA_COPY)

if (UNIX)
  set (HAVE_PROC 1)
else (UNIX)
  set (HAVE_PROC)
endif (UNIX)

# dependencies

# parallel programming
# enabling OpenMP is supposedly enough to make the compiler link with
# the appropriate libraries
find_package (OpenMP ${ERT_QUIET})
list (APPEND ERT_LIBRARIES ${OpenMP_LIBRARIES})
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

# threading library (search for this *after* OpenMP
set (CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package (Threads ${ERT_QUIET})
if (CMAKE_USE_PTHREADS_INIT)
  list (APPEND ERT_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif (CMAKE_USE_PTHREADS_INIT)

# compression library
find_package (ZLIB ${ERT_QUIET})
if (ZLIB_FOUND)
  list (APPEND ERT_INCLUDE_DIRS ${ZLIB_INCLUDE_DIRS})
  list (APPEND ERT_LIBRARIES ${ZLIB_LIBRARIES})
endif (ZLIB_FOUND)

# numerics
find_package (BLAS ${ERT_QUIET})
if (BLAS_FOUND)
  list (APPEND ERT_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
  list (APPEND ERT_LIBRARIES ${BLAS_LIBRARIES})
  list (APPEND ERT_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
endif (BLAS_FOUND)
find_package (LAPACK ${ERT_QUIET})
if (LAPACK_FOUND)
  list (APPEND ERT_INCLUDE_DIRS ${LAPACK_INCLUDE_DIRS})
  list (APPEND ERT_LIBRARIES ${LAPACK_LIBRARIES})
  list (APPEND ERT_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
endif (LAPACK_FOUND)

# math library (should exist on all unices; automatically linked on Windows)
if (UNIX)
  find_library (MATH_LIBRARY
	NAMES "m"
	)
  list (APPEND ERT_LIBRARIES ${MATH_LIBRARY})
endif (UNIX)

# since OpenMP often implies pthreads, we need to tidy up
# (last instance of library must be left standing, thus reversing that
# list before removing duplicates)
list (REMOVE_DUPLICATES ERT_INCLUDE_DIRS)
list (REVERSE ERT_LIBRARIES)
list (REMOVE_DUPLICATES ERT_LIBRARIES)
list (REVERSE ERT_LIBRARIES)

# linker flags may not be specified at all
if (DEFINED ERT_LINKER_FLAGS)
  list (REMOVE_DUPLICATES ERT_LINKER_FLAGS)
endif (DEFINED ERT_LINKER_FLAGS)

# see if we can compile a minimum example
# CMake logical test doesn't handle lists (sic)
if (NOT ERT_LIBRARIES MATCHES "-NOTFOUND")
  cmake_push_check_state ()
  include (CheckCSourceCompiles)
  set (CMAKE_REQUIRED_INCLUDES ${ERT_INCLUDE_DIR})
  set (CMAKE_REQUIRED_LIBRARIES ${ERT_LIBRARIES})
  check_c_source_compiles (
	"#include <ecl_util.h>
int main (void) {
  int sz;
  sz = ecl_util_get_sizeof_ctype (ECL_INT_TYPE);
  return 0;
}" HAVE_ERT)
  cmake_pop_check_state ()
else (NOT ERT_LIBRARIES MATCHES "-NOTFOUND")
  set (HAVE_ERT 0)
endif (NOT ERT_LIBRARIES MATCHES "-NOTFOUND")

# if the test program didn't compile, but was required to do so, bail
# out now and display an error; otherwise limp on
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (ERT
  DEFAULT_MSG
  ERT_INCLUDE_DIR ERT_LIBRARY HAVE_ERT
  )
