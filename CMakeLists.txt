cmake_minimum_required (VERSION 2.6)
project (opm-core)

enable_language(Fortran)



set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.36.0 COMPONENTS date_time filesystem system unit_test_framework)
find_package(UmfPack)
find_package(LAPACK)

include_directories("${PROJECT_SOURCE_DIR} ${Boost_INCLUDE_DIRS}")


# The opmcore library
add_library(opmcore
opm/core/eclipse/EclipseGridInspector.cpp
opm/core/eclipse/EclipseGridParser.cpp
opm/core/fluid/blackoil/BlackoilPvtProperties.cpp
opm/core/fluid/blackoil/SinglePvtDead.cpp
opm/core/fluid/blackoil/SinglePvtLiveGas.cpp
opm/core/fluid/blackoil/SinglePvtLiveOil.cpp
opm/core/fluid/blackoil/SinglePvtInterface.cpp
opm/core/fluid/BlackoilPropertiesBasic.cpp
opm/core/fluid/BlackoilPropertiesFromDeck.cpp
opm/core/fluid/IncompPropertiesBasic.cpp
opm/core/fluid/IncompPropertiesFromDeck.cpp
opm/core/fluid/PvtPropertiesBasic.cpp
opm/core/fluid/PvtPropertiesIncompFromDeck.cpp
opm/core/fluid/RockBasic.cpp
opm/core/fluid/RockFromDeck.cpp
opm/core/fluid/SaturationPropsBasic.cpp
opm/core/fluid/SaturationPropsFromDeck.cpp
opm/core/utility/MonotCubicInterpolator.cpp
opm/core/utility/parameters/Parameter.cpp
opm/core/utility/parameters/ParameterGroup.cpp
opm/core/utility/parameters/ParameterTools.cpp
opm/core/utility/parameters/ParameterXML.cpp
opm/core/utility/parameters/tinyxml/tinystr.cpp
opm/core/utility/parameters/tinyxml/tinyxml.cpp
opm/core/utility/parameters/tinyxml/tinyxmlerror.cpp
opm/core/utility/parameters/tinyxml/tinyxmlparser.cpp
opm/core/utility/cart_grid.c
opm/core/utility/cpgpreprocess/geometry.c
opm/core/utility/cpgpreprocess/preprocess.c
opm/core/utility/cpgpreprocess/readvector.cpp
opm/core/utility/cpgpreprocess/cgridinterface.c
opm/core/utility/cpgpreprocess/sparsetable.c
opm/core/utility/cpgpreprocess/facetopology.c
opm/core/utility/cpgpreprocess/uniquepoints.c
opm/core/utility/StopWatch.cpp
opm/core/utility/newwells.c
opm/core/utility/writeVtkData.cpp
opm/core/GridManager.cpp
opm/core/linalg/sparse_sys.c
opm/core/pressure/cfsh.c
opm/core/pressure/flow_bc.c
opm/core/pressure/well.c
opm/core/pressure/fsh_common_impl.c
opm/core/pressure/fsh.c
opm/core/pressure/tpfa/ifs_tpfa.c
opm/core/pressure/tpfa/compr_source.c
opm/core/pressure/tpfa/cfs_tpfa.c
opm/core/pressure/tpfa/compr_bc.c
opm/core/pressure/tpfa/compr_quant.c
opm/core/pressure/tpfa/compr_quant_general.c
opm/core/pressure/tpfa/cfs_tpfa_residual.c
opm/core/pressure/tpfa/trans_tpfa.c
opm/core/pressure/msmfem/coarse_conn.c
opm/core/pressure/msmfem/partition.c
opm/core/pressure/msmfem/hash_set.c
opm/core/pressure/msmfem/ifsh_ms.c
opm/core/pressure/msmfem/dfs.c
opm/core/pressure/msmfem/coarse_sys.c
opm/core/pressure/ifsh.c
opm/core/pressure/IncompTpfa.cpp
opm/core/pressure/mimetic/mimetic.c
opm/core/pressure/mimetic/hybsys_global.c
opm/core/pressure/mimetic/hybsys.c
opm/core/transport/spu_explicit.c
opm/core/transport/spu_implicit.c
opm/core/transport/transport_source.c
opm/core/transport/reorder/TransportModelInterface.cpp
opm/core/transport/reorder/TransportModelTwophase.cpp
opm/core/transport/reorder/reordersequence.cpp
opm/core/transport/reorder/nlsolvers.c
opm/core/transport/reorder/tarjan.c
opm/core/linalg/call_umfpack.c
opm/core/pressure/IncompTpfa.cpp
)

target_link_libraries(opmcore
  ${UMFPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
 -lcholmod -lcamd -lccolamd -lmetis
)

add_executable(spu_2p examples/spu_2p.cpp)

target_link_libraries(spu_2p
  opmcore
)

set_target_properties(opmcore spu_2p PROPERTIES COMPILE_FLAGS -m64 LINKER_LANGUAGE CXX LINK_FLAGS -m64 )
