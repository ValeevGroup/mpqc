cmake_minimum_required(VERSION 2.8.11)
project(HeirChem)

set(HeirChem_VERSION_MAJOR 0)
set(HeirChem_VERSION_MINOR 0)

option(RELEASE ON)

set(CXX_FLAGS "-std=c++11")
if(RELEASE)
    set(CXX_FLAGS "${CXX_FLAGS} -DNDEBUG -mavx -O3")
    message(STATUS "Performing release build ${CXX_FLAGS}")
ELSEIF(DEBUG_RELEASE)
    set(CXX_FLAGS "${CXX_FLAGS} -O2 -g")
    message(STATUS "Performing debug-release build ${CXX_FLAGS}")
ELSE()
    set(CXX_FLAGS "${CXX_FLAGS} -g -Wall -Wextra -Werror")
    message(STATUS "Performing debug build ${CXX_FLAGS}")
ENDIF()

set(CMAKE_CXX_FLAGS ${CXX_FLAGS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")

FIND_PACKAGE(MPI REQUIRED)
SET(EXTRA_FLAGS "${MPI_CXX_COMPILE_FLAGS}")

FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(MPQC REQUIRED)
FIND_PACKAGE(TILEDARRAY REQUIRED)

INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${MPQC_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${TILEDARRAY_INCLUDE_DIRS})
LIST(APPEND HeirChem_LIBRARIES ${LAPACK_LIBRARIES})
LIST(APPEND HeirChem_LIBRARIES ${MPI_CXX_LIBRARIES})
LIST(APPEND HeirChem_LIBRARIES ${TBB_LIBRARIES})
LIST(APPEND HeirChem_LIBRARIES ${MPQC_LIBRARIES})
LIST(APPEND HeirChem_LIBRARIES ${TILEDARRAY_LIBRARIES})




INCLUDE_DIRECTORIES(${HeirChem_SOURCE_DIR}/molecule)
INCLUDE_DIRECTORIES(${HeirChem_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${HeirChem_SOURCE_DIR}/matrix)

add_subdirectory(${HeirChem_SOURCE_DIR}/molecule)
add_subdirectory(${HeirChem_SOURCE_DIR}/include)
add_subdirectory(${HeirChem_SOURCE_DIR}/matrix)

LINK_DIRECTORIES(${HeirChem_BINARY_DIR}/molecule)
LIST(APPEND HeirChem_LIBRARIES molecule)

message(STATUS "Using libraries ${HeirChem_LIBRARIES}")

set(HeirChem_Sources "main.cpp")

add_executable(${PROJECT_NAME} ${HeirChem_Sources}) 
target_link_libraries(${PROJECT_NAME}
  ${HeirChem_LIBRARIES}
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
