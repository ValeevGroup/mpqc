cmake_minimum_required(VERSION 3.0.0)
project(HeirChem)

set(HeirChem_VERSION_MAJOR 0)
set(HeirChem_VERSION_MINOR 0)

option(RELEASE ON)

if(RELEASE)
    set(CXX_FLAGS "-O3 -std=c++11" CACHE STRING "CXX flags")
ELSEIF(DEBUG_RELEASE)
    set(CXX_FLAGS "-O3 -g -Rpass-analysis=loop-vectorize -std=c++11" CACHE STRING "CXX flags")
else()
    set(CXX_FLAGS "-g -Wall -Wextra -std=c++11" CACHE STRING "CXX flags")
endif()


# Now append the extra flags
set(CMAKE_CXX_FLAGS "${CXX_FLAGS} ${EXTRA_FLAGS}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")
FIND_PACKAGE(MPI REQUIRED)
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(dmhm)
FIND_PACKAGE(MPQC)

if(NOT MPI_CXX_FOUND)
  message(FATAL_ERROR "MPI C++ compiler was not found and is required")
endif()
include_directories(${MPI_CXX_INCLUDE_PATH})

IF(${LAPACK_FOUND})
  message(STATUS "LAPACK FOUND")
  INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIRS})
  list(APPEND HeirChem_LIBRARIES ${LAPACK_LIBRARIES})
ENDIF()

IF(${EIGEN_FOUND})
  message(STATUS "EIGEN FOUND")
  INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
ENDIF()

if(${DMHM_FOUND})
  message(STATUS "DMHM FOUND")
  include_directories( ${DMHM_INCLUDE_DIRS} )
  list(APPEND HeirChem_LIBRARIES ${DMHM_LIBRARIES})
endif()

if(${TBB_FOUND})
  include_directories( ${TBB_INCLUDE_DIRS} )
  list(APPEND HeirChem_LIBRARIES ${TBB_LIBRARIES})
endif()

if(${MPQC_FOUND})
  include_directories( ${MPQC_INCLUDE_DIRS} )
  list(APPEND HeirChem_LIBRARIES ${MPQC_LIBRARIES})
endif()

set(EXTRA_FLAGS "${MPI_CXX_COMPILE_FLAGS}")


INCLUDE_DIRECTORIES(${HeirChem_SOURCE_DIR}/molecule)
INCLUDE_DIRECTORIES(${HeirChem_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${HeirChem_SOURCE_DIR}/matrix)

add_subdirectory(${HeirChem_SOURCE_DIR}/molecule)
add_subdirectory(${HeirChem_SOURCE_DIR}/include)

LINK_DIRECTORIES(${HeirChem_BINARY_DIR}/molecule)

LIST(APPEND HeirChem_LIBRARIES molecule)

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
