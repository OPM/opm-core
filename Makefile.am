ACLOCAL_AMFLAGS = -I m4

# Recurse to build examples after libraries
SUBDIRS = . tests examples

# Declare libraries
lib_LTLIBRARIES = libopmcore.la

libopmcore_la_SOURCES = 				\
opm/core/eclipse/EclipseGridInspector.cpp               \
opm/core/eclipse/EclipseGridParser.cpp                  \
opm/core/fluid/blackoil/BlackoilPVT.cpp                 \
opm/core/fluid/blackoil/MiscibilityDead.cpp             \
opm/core/fluid/blackoil/MiscibilityLiveGas.cpp          \
opm/core/fluid/blackoil/MiscibilityLiveOil.cpp          \
opm/core/fluid/blackoil/MiscibilityProps.cpp            \
opm/core/utility/MonotCubicInterpolator.cpp             \
opm/core/utility/parameters/Parameter.cpp               \
opm/core/utility/parameters/ParameterGroup.cpp          \
opm/core/utility/parameters/ParameterTools.cpp          \
opm/core/utility/parameters/ParameterXML.cpp            \
opm/core/utility/parameters/tinyxml/tinystr.cpp         \
opm/core/utility/parameters/tinyxml/tinyxml.cpp         \
opm/core/utility/parameters/tinyxml/tinyxmlerror.cpp    \
opm/core/utility/parameters/tinyxml/tinyxmlparser.cpp   \
opm/core/utility/cart_grid.c				\
opm/core/utility/cpgpreprocess/geometry.c		\
opm/core/utility/cpgpreprocess/preprocess.c		\
opm/core/utility/cpgpreprocess/readvector.cpp		\
opm/core/utility/cpgpreprocess/cgridinterface.c		\
opm/core/utility/cpgpreprocess/sparsetable.c		\
opm/core/utility/cpgpreprocess/facetopology.c		\
opm/core/utility/cpgpreprocess/uniquepoints.c		\
opm/core/utility/StopWatch.cpp                          \
opm/core/linalg/sparse_sys.c				\
opm/core/pressure/cfsh.c				\
opm/core/pressure/flow_bc.c				\
opm/core/pressure/well.c				\
opm/core/pressure/fsh_common_impl.c			\
opm/core/pressure/fsh.c			                \
opm/core/pressure/tpfa/ifs_tpfa.c                       \
opm/core/pressure/tpfa/compr_source.c                   \
opm/core/pressure/tpfa/cfs_tpfa.c                       \
opm/core/pressure/tpfa/compr_bc.c                       \
opm/core/pressure/tpfa/compr_quant.c                    \
opm/core/pressure/tpfa/compr_quant_general.c            \
opm/core/pressure/tpfa/cfs_tpfa_residual.c              \
opm/core/pressure/tpfa/trans_tpfa.c                     \
opm/core/pressure/msmfem/coarse_conn.c                  \
opm/core/pressure/msmfem/partition.c                    \
opm/core/pressure/msmfem/hash_set.c                     \
opm/core/pressure/msmfem/ifsh_ms.c                      \
opm/core/pressure/msmfem/dfs.c                          \
opm/core/pressure/msmfem/coarse_sys.c                   \
opm/core/pressure/ifsh.c                                \
opm/core/pressure/mimetic/mimetic.c                     \
opm/core/pressure/mimetic/hybsys_global.c               \
opm/core/pressure/mimetic/hybsys.c                      \
opm/core/transport/spu_explicit.c                       \
opm/core/transport/spu_implicit.c                       \
opm/core/transport/transport_source.c

nobase_include_HEADERS =                                \
opm/core/eclipse/CornerpointChopper.hpp                 \
opm/core/eclipse/EclipseGridInspector.hpp               \
opm/core/eclipse/EclipseGridParser.hpp                  \
opm/core/eclipse/EclipseGridParserHelpers.hpp           \
opm/core/eclipse/EclipseUnits.hpp                       \
opm/core/eclipse/SpecialEclipseFields.hpp               \
opm/core/fluid/BlackoilFluid.hpp                        \
opm/core/fluid/SimpleFluid2p.hpp                        \
opm/core/fluid/blackoil/BlackoilDefs.hpp                \
opm/core/fluid/blackoil/BlackoilPVT.hpp                 \
opm/core/fluid/blackoil/FluidMatrixInteractionBlackoil.hpp \
opm/core/fluid/blackoil/FluidStateBlackoil.hpp          \
opm/core/fluid/blackoil/MiscibilityDead.hpp             \
opm/core/fluid/blackoil/MiscibilityLiveGas.hpp          \
opm/core/fluid/blackoil/MiscibilityLiveOil.hpp          \
opm/core/fluid/blackoil/MiscibilityProps.hpp            \
opm/core/fluid/blackoil/MiscibilityWater.hpp            \
opm/core/utility/Average.hpp                            \
opm/core/utility/ErrorMacros.hpp                        \
opm/core/utility/Factory.hpp                            \
opm/core/utility/MonotCubicInterpolator.hpp             \
opm/core/utility/parameters/Parameter.hpp               \
opm/core/utility/parameters/ParameterGroup.hpp          \
opm/core/utility/parameters/ParameterGroup_impl.hpp     \
opm/core/utility/parameters/ParameterMapItem.hpp        \
opm/core/utility/parameters/ParameterRequirement.hpp    \
opm/core/utility/parameters/ParameterStrings.hpp        \
opm/core/utility/parameters/ParameterTools.hpp          \
opm/core/utility/parameters/ParameterXML.hpp            \
opm/core/utility/parameters/tinyxml/tinystr.h           \
opm/core/utility/parameters/tinyxml/tinyxml.h           \
opm/core/utility/RootFinders.hpp                        \
opm/core/utility/SparseTable.hpp                        \
opm/core/utility/SparseVector.hpp                       \
opm/core/utility/StopWatch.hpp                          \
opm/core/utility/UniformTableLinear.hpp                 \
opm/core/utility/Units.hpp                              \
opm/core/utility/buildUniformMonotoneTable.hpp          \
opm/core/utility/linInt.hpp                             \
opm/core/utility/linearInterpolation.hpp                \
opm/core/utility/cpgpreprocess/readvector.hpp           \
opm/core/utility/cpgpreprocess/uniquepoints.h           \
opm/core/utility/cpgpreprocess/preprocess.h             \
opm/core/utility/cpgpreprocess/sparsetable.h            \
opm/core/utility/cpgpreprocess/cgridinterface.h         \
opm/core/utility/cpgpreprocess/geometry.h               \
opm/core/utility/cpgpreprocess/facetopology.h           \
opm/core/utility/cpgpreprocess/grid.h                   \
opm/core/utility/cpgpreprocess/grdecl.h                 \
opm/core/utility/cart_grid.h                            \
opm/core/well.h                                         \
opm/core/linalg/sparse_sys.h                            \
opm/core/linalg/blas_lapack.h                           \
opm/core/grid.h                                         \
opm/core/GridAdapter.hpp                                \
opm/core/pressure/fsh.h                                 \
opm/core/pressure/HybridPressureSolver.hpp              \
opm/core/pressure/TPFAPressureSolver.hpp                \
opm/core/pressure/tpfa/compr_quant.h                    \
opm/core/pressure/tpfa/compr_source.h                   \
opm/core/pressure/tpfa/compr_bc.h                       \
opm/core/pressure/tpfa/cfs_tpfa.h                       \
opm/core/pressure/tpfa/compr_quant_general.h            \
opm/core/pressure/tpfa/trans_tpfa.h                     \
opm/core/pressure/tpfa/cfs_tpfa_residual.h              \
opm/core/pressure/tpfa/ifs_tpfa.h                       \
opm/core/pressure/TPFACompressiblePressureSolver.hpp    \
opm/core/pressure/msmfem/coarse_sys.h                   \
opm/core/pressure/msmfem/hash_set.h                     \
opm/core/pressure/msmfem/ifsh_ms.h                      \
opm/core/pressure/msmfem/dfs.h                          \
opm/core/pressure/msmfem/partition.h                    \
opm/core/pressure/msmfem/coarse_conn.h                  \
opm/core/pressure/flow_bc.h                             \
opm/core/pressure/fsh_common_impl.h                     \
opm/core/pressure/mimetic/hybsys_global.h               \
opm/core/pressure/mimetic/hybsys.h                      \
opm/core/pressure/mimetic/mimetic.h                     \
opm/core/transport/ImplicitTransport.hpp                \
opm/core/transport/CSRMatrixUmfpackSolver.hpp           \
opm/core/transport/CSRMatrixBlockAssembler.hpp          \
opm/core/transport/SimpleFluid2pWrapper.hpp             \
opm/core/transport/SinglePointUpwindTwoPhase.hpp        \
opm/core/transport/NormSupport.hpp                      \
opm/core/transport/spu_explicit.h                       \
opm/core/transport/JacobianSystem.hpp                   \
opm/core/transport/spu_implicit.h                       \
opm/core/transport/ImplicitAssembly.hpp                 \
opm/core/transport/transport_source.h

if UMFPACK
libopmcore_la_SOURCES  += opm/core/linalg/call_umfpack.c
nobase_include_HEADERS += opm/core/linalg/call_umfpack.h
endif
