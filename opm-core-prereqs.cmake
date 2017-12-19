# defines that must be present in config.h for our headers
set (opm-core_CONFIG_VAR
  HAVE_ERT
  HAVE_SUITESPARSE_UMFPACK_H
  HAVE_DUNE_ISTL
  HAVE_MPI
  HAVE_PETSC
  DUNE_ISTL_VERSION_MAJOR
  DUNE_ISTL_VERSION_MINOR
  DUNE_ISTL_VERSION_REVISION
  )

# dependencies
set (opm-core_DEPS
  # compile with C99 support if available
  "C99"
  # compile with C++0x/11 support if available
  "CXX11Features REQUIRED"
  # various runtime library enhancements
  "Boost 1.44.0
    COMPONENTS date_time filesystem system unit_test_framework REQUIRED"
  # matrix library
  "BLAS REQUIRED"
  "LAPACK REQUIRED"
  # Tim Davis' SuiteSparse archive
  "SuiteSparse COMPONENTS umfpack"
  # solver
  "SuperLU"
	# Eclipse I/O tools
  "ecl REQUIRED"
  # Look for MPI support
  "MPI"
  # PETSc numerical backend
  "PETSc"
  # DUNE dependency
  "dune-common"
  "dune-istl"
  "opm-common REQUIRED"
  # Parser library for ECL-type simulation models
  "opm-parser REQUIRED"
  # the code which implements the material laws
  "opm-material REQUIRED"
  # the code which implements the output routines
  "opm-output REQUIRED"
  # the code which implements grids
  "opm-grid REQUIRED"
  )

find_package_deps(opm-core)
