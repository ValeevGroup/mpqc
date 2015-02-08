cmake_minimum_required(VERSION 2.8.11)
project(TileClusterChem)

set(TileClusterChem_VERSION_MAJOR 0)
set(TileClusterChem_VERSION_MINOR 1)
SET(CMAKE_COLOR_MAKEFILE ON)

option(RELEASE OFF)
option(DEBUG_RELEASE OFF)
enable_testing()
option(PERFORM_UNIT_TESTS ON)

set(CXX_FLAGS "${CXX_FLAGS} -std=c++11")
if(RELEASE)
    set(CXX_FLAGS "${CXX_FLAGS} -DNDEBUG -mavx -O3")
    message(STATUS "Performing release build ${CXX_FLAGS}")
ELSEIF(DEBUG_RELEASE)
    set(CXX_FLAGS "${CXX_FLAGS} -O2 -mavx -g ")
    message(STATUS "Performing debug-release build ${CXX_FLAGS}")
ELSEIF(FAST_BUILD)
    set(CXX_FLAGS "${CXX_FLAGS} -O0 -g -Wall -Wextra -Wno-unused-parameter")
    message(STATUS "Performing debug-release build ${CXX_FLAGS}")
ELSE()
    set(CXX_FLAGS "${CXX_FLAGS} -g -Wall -Wextra -Wno-unused-parameter")
    message(STATUS "Performing debug build ${CXX_FLAGS}")
ENDIF()

set(CMAKE_CXX_FLAGS ${CXX_FLAGS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/cmake/modules/")

FIND_PACKAGE(MPI REQUIRED)
SET(EXTRA_FLAGS "${MPI_CXX_COMPILE_FLAGS}")

include(external/lapack.cmake)

FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(LIBINT REQUIRED)
FIND_PACKAGE(BTAS REQUIRED)
FIND_PACKAGE(TILEDARRAY REQUIRED)

INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${TILEDARRAY_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LIBINT_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${BTAS_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIRS})
LIST(APPEND TileClusterChem_LIBRARIES ${LAPACK_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${MPI_CXX_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${TBB_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${LIBINT_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${TILEDARRAY_LIBRARIES})

add_subdirectory(${TileClusterChem_SOURCE_DIR}/include)

add_subdirectory(${TileClusterChem_SOURCE_DIR}/molecule)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/molecule)
LIST(APPEND TileClusterChem_LIBRARIES molecule)

add_subdirectory(${TileClusterChem_SOURCE_DIR}/tensor)

add_subdirectory(${TileClusterChem_SOURCE_DIR}/basis)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/basis)
LIST(APPEND TileClusterChem_LIBRARIES basis)

add_subdirectory(${TileClusterChem_SOURCE_DIR}/integrals)

message(STATUS "Using libraries ${TileClusterChem_LIBRARIES}")

set(TileClusterChem_Sources "main.cpp")

IF(PERFORM_UNIT_TESTS)
  add_subdirectory(${TileClusterChem_SOURCE_DIR}/tests)
ENDIF(PERFORM_UNIT_TESTS)

add_executable(${PROJECT_NAME} ${TileClusterChem_Sources}) 
target_link_libraries(${PROJECT_NAME}
  ${TileClusterChem_LIBRARIES}
)

add_executable(sparse_integral_example EXCLUDE_FROM_ALL sparse_integral_example.cpp) 
target_link_libraries(sparse_integral_example
  ${TileClusterChem_LIBRARIES}
)

add_executable(integral_space_savings EXCLUDE_FROM_ALL integral_space_savings.cpp) 
target_link_libraries(integral_space_savings
  ${TileClusterChem_LIBRARIES}
)

add_executable(integral_correctness_test EXCLUDE_FROM_ALL integral_correctness_test.cpp) 
target_link_libraries(integral_correctness_test
  ${TileClusterChem_LIBRARIES}
)

add_executable(sqrt_inverse_example EXCLUDE_FROM_ALL sqrt_inverse_example.cpp) 
target_link_libraries(sqrt_inverse_example
  ${TileClusterChem_LIBRARIES}
)

add_executable(eri3_low_rank_example EXCLUDE_FROM_ALL eri3_low_rank_example.cpp) 
target_link_libraries(eri3_low_rank_example
  ${TileClusterChem_LIBRARIES}
)

# Makes QT creator show the following file types
file(GLOB_RECURSE QT_CREATOR_SRC
  "*.h"
  "*.hpp"
  "*.cpp"
  "*.c"
  "*.txt"
  "*.in"
  "*.cc"
  "TODO"
  ".clang-format"
  "README"
  )

add_library(qt_creator_get_sources EXCLUDE_FROM_ALL ${QT_CREATOR_SRC})

set_target_properties(qt_creator_get_sources PROPERTIES LINKER_LANGUAGE CXX)
