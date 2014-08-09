project(HeirChem)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")
FIND_PACKAGE(Eigen3)

add_library(Eigen3 INTERFACE)
target_include_directories(Eigen3 INTERFACE
  ${EIGEN_INCLUCE_DIRS}
)

set(HeirChem_Sources "main.cpp")

set(directories
  molecule
  include
)

include_directories(
  ${directories}
  ${EIGEN_INCLUDE_DIRS}
  )

add_subdirectory(${directories})

message(STATUS "${HeirChem_Sources}")

add_executable(${PROJECT_NAME} ${HeirChem_Sources}) 
target_link_libraries(${PROJECT_NAME} Eigen3)

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
  )

add_library(qt_creator_get_sources EXCLUDE_FROM_ALL ${QT_CREATOR_SRC})

set_target_properties(qt_creator_get_sources PROPERTIES LINKER_LANGUAGE CXX)
