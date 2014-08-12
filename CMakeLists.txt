cmake_minimum_required(VERSION 3.0.0)
project(HeirChem)

set(HeirChem_VERSION_MAJOR 0)
set(HeirChem_VERSION_MINOR 0)

option(RELEASE ON)

set(CMAKE_CXX_COMPILER "clang++")

if(RELEASE)
    set(CXX_FLAGS "-O3 -std=c++11" CACHE STRING "CXX flags")
else()
    set(CXX_FLAGS "-g -Wall -Wextra -std=c++11" CACHE STRING "CXX flags")
endif()

find_package(MPI)
if(NOT MPI_CXX_FOUND)
  message(FATAL_ERROR "MPI C++ compiler was not found and is required")
endif()
include_directories(${MPI_CXX_INCLUDE_PATH})
set(EXTRA_FLAGS "${MPI_CXX_COMPILE_FLAGS}")

# Now append the extra flags
set(CMAKE_CXX_FLAGS "${CXX_FLAGS} ${EXTRA_FLAGS}")

if(NOT MATH_LIBS)
  set(MATH_DESC "BLAS/LAPACK link flags")
  # Look for default BLAS/LAPACK
  set(REFERENCE_REQUIRED LAPACK BLAS)
  find_library(BLAS_LIB
               NAMES blas blas.a blas.lib
               PATHS ${REFERENCE_ROOT})
  find_library(LAPACK_LIB
               NAMES lapack lapack.a lapack.lib
                     reflapack reflapack.a reflapack.lib
               PATHS ${REFERENCE_ROOT})
  set(REFERENCE_FOUND ON)
  set(MATH_LIBS "")
  foreach(NAME ${REFERENCE_REQUIRED})
    if( ${NAME}_LIB )
      message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
      list(APPEND MATH_LIBS ${${NAME}_LIB})
    else()
      message(STATUS "Could not find ${NAME}_LIB")
      set(REFERENCE_FOUND OFF)
    endif()
  endforeach()
  message(STATUS "REFERENCE_FOUND=${REFERENCE_FOUND}")
  if(REFERENCE_FOUND)
    message(STATUS "WARNING: Using reference BLAS/LAPACK.")
    message(STATUS "MATH_LIBS=${MATH_LIBS}")
  else()
    set(MATH_LIBS "" CACHE STRING ${MATH_DESC})
    message(FATAL_ERROR
      "Could not find BLAS/LAPACK libs. Please provide the root directory of MKL with -DMKL_ROOT, the directory of reference implementations with -DREFERENCE_ROOT, or manually specify all math libraries with -DMATH_LIBS. There are numerous idiosyncratic library dependencies for BLAS/LAPACK, so you will almost certainly need to manually specify -DMATH_LIBS.")
  endif()
  # Append the standard math libraries to the link list.
  list(APPEND MATH_LIBS m)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")
FIND_PACKAGE(Eigen3)
FIND_PACKAGE(TBB)
FIND_PACKAGE(dmhm)
FIND_PACKAGE(MPQC)

IF(${EIGEN_FOUND})
  message(STATUS "EIGEN FOUND")
  INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
ENDIF()

if(${DMHM_FOUND})
  include_directories(
    ${DMHM_INCLUDE_DIRS}
    )
  list(APPEND HeirChem_LIBRARIES ${DMHM_PATH}/lib/lib${DMHM_LIBRARIES}.a)
endif()

if(${TBB_FOUND})
  include_directories(
    ${TBB_INCLUDE_DIRS}
    )
  list(APPEND HeirChem_LIBRARIES ${TBB_LIBRARIES})
endif()

if(${MPQC_FOUND})
  message(STATUS "${MPQC_FOUND}")
  message(STATUS "${MPQC_INCLUDE_DIRS}")
  include_directories(
    ${MPQC_INCLUDE_DIRS}
    )
  list(APPEND HeirChem_LIBRARIES ${MPQC_LIBRARIES})
endif()


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
