cmake_minimum_required (VERSION 2.6)
project (opm-core)

enable_language(Fortran)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})


set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.36.0 COMPONENTS date_time filesystem system unit_test_framework)
find_package(UmfPack)
find_package(LAPACK)

include_directories(${PROJECT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})


# The opmcore library
FILE(GLOB_RECURSE C_FILES_CORE "opm/core/*.c")
FILE(GLOB_RECURSE CPP_FILES_CORE "opm/core/*.cpp")
FILE(GLOB_RECURSE REMOVE_FILES   "processgrid.c")
FILE(GLOB_RECURSE REMOVE_FILESMX   "mx*.c")
FILE(GLOB_RECURSE REMOVE_FILESAGMG   "*AGMG.cpp" "*test*")
list(REMOVE_ITEM C_FILES_CORE ${REMOVE_FILES} ${REMOVE_FILESMX} )
list(REMOVE_ITEM CPP_FILES_CORE ${REMOVE_FILES} ${REMOVE_FILESAGMG})
add_library(opmcore ${C_FILES_CORE} ${CPP_FILES_CORE} )


target_link_libraries(opmcore
  ${UMFPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES} ${Boost_LIBRARIES}
 -lcholmod -lcamd -lccolamd -lmetis -ldunecommon
)

FILE(GLOB CPP_EXAMPLES "examples/*.cpp")
FILE(GLOB CPP_tests "tests/*.cpp")
add_executable(spu_2p examples/spu_2p.cpp)
#add_executable(sim_2p_incomp_reorder examples/sim_2p_incomp_reorder.cpp)
#add_executable(sim_water_oil examples/sim_water_oil.cpp)

add_executable(pvt_test tests/pvt_test.cpp)
add_executable(relperm_test tests/relperm_test.cpp)

target_link_libraries(spu_2p
  opmcore
)
target_link_libraries(pvt_test
  opmcore
)
target_link_libraries(relperm_test
  opmcore
)

set_target_properties(opmcore spu_2p PROPERTIES COMPILE_FLAGS -m64 LINKER_LANGUAGE CXX LINK_FLAGS -m64)
set_target_properties(opmcore pvt_test PROPERTIES COMPILE_FLAGS -m64 LINKER_LANGUAGE CXX LINK_FLAGS -m64)
set_target_properties(opmcore relperm_test PROPERTIES COMPILE_FLAGS -m64 LINKER_LANGUAGE CXX LINK_FLAGS -m64)
