cmake_minimum_required(VERSION 2.8.11)
project(TileClusterChem)

set(TileClusterChem_VERSION_MAJOR 0)
set(TileClusterChem_VERSION_MINOR 1)

option(RELEASE OFF)
option(DEBUG_RELEASE OFF)
option(PERFORM_UNIT_TESTS ON)

set(CXX_FLAGS "-std=c++1y")
if(RELEASE)
    set(CXX_FLAGS "${CXX_FLAGS} -DNDEBUG -mavx -O3")
    message(STATUS "Performing release build ${CXX_FLAGS}")
ELSEIF(DEBUG_RELEASE)
    set(CXX_FLAGS "${CXX_FLAGS} -O2 -mavx -g")
    message(STATUS "Performing debug-release build ${CXX_FLAGS}")
ELSE()
    set(CXX_FLAGS "${CXX_FLAGS} -fno-inline -g -Wall -Wextra -Wno-unused-parameter")
    message(STATUS "Performing debug build ${CXX_FLAGS}")
ENDIF()

set(CMAKE_CXX_FLAGS ${CXX_FLAGS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/cmake/modules/")

FIND_PACKAGE(MPI REQUIRED)
SET(EXTRA_FLAGS "${MPI_CXX_COMPILE_FLAGS}")

if(USER_MATH_LIBS)
  LIST(APPEND EXTRA_FLAGS ${USER_MATH_CXX_FLAGS})
  Set(LAPACK_INCLUDE_DIRS ${USER_MATH_INCLUDE_DIRS})
  Set(LAPACK_LIBRARIES ${USER_MATH_LIBS})
endif(USER_MATH_LIBS)

if(NOT USER_MATH_LIBS)
  FIND_PACKAGE(LAPACK REQUIRED)
endif(NOT USER_MATH_LIBS)

FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(MPQC)
FIND_PACKAGE(LIBINT REQUIRED)
FIND_PACKAGE(TILEDARRAY REQUIRED)

INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${MPQC_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${TILEDARRAY_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LIBINT_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIRS})
LIST(APPEND TileClusterChem_LIBRARIES ${LAPACK_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${MPI_CXX_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${TBB_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${MPQC_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${LIBINT_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${TILEDARRAY_LIBRARIES})


INCLUDE_DIRECTORIES(${TileClusterChem_SOURCE_DIR}/molecule)
INCLUDE_DIRECTORIES(${TileClusterChem_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${TileClusterChem_SOURCE_DIR}/matrix)


add_subdirectory(${TileClusterChem_SOURCE_DIR}/include)
add_subdirectory(${TileClusterChem_SOURCE_DIR}/molecule)
add_subdirectory(${TileClusterChem_SOURCE_DIR}/matrix)
add_subdirectory(${TileClusterChem_SOURCE_DIR}/basis)
add_subdirectory(${TileClusterChem_SOURCE_DIR}/integrals)

IF(PERFORM_UNIT_TESTS)
  add_subdirectory(${TileClusterChem_SOURCE_DIR}/tests)
ENDIF(PERFORM_UNIT_TESTS)

LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/molecule)
LIST(APPEND TileClusterChem_LIBRARIES molecule)

message(STATUS "Using libraries ${TileClusterChem_LIBRARIES}")

set(TileClusterChem_Sources "main.cpp")

add_executable(${PROJECT_NAME} ${TileClusterChem_Sources}) 
target_link_libraries(${PROJECT_NAME}
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
