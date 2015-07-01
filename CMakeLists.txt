cmake_minimum_required(VERSION 2.8.11)
project(TileClusterChem)

set(TileClusterChem_VERSION_MAJOR 0)
set(TileClusterChem_VERSION_MINOR 1)
SET(CMAKE_COLOR_MAKEFILE ON)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/cmake/modules/")

FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(Boost REQUIRED)
FIND_PACKAGE(LIBINT REQUIRED)
FIND_PACKAGE(BTAS REQUIRED)
FIND_PACKAGE(TiledArray CONFIG HINTS ${TiledArray_CONFIG_DIR} REQUIRED)

INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${TiledArray_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LIBINT_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${BTAS_INCLUDE_DIRS})

LIST(APPEND CMAKE_CXX_FLAGS ${TiledArray_COMPILE_FLAGS})

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

LIST(APPEND TileClusterChem_LIBRARIES ${TBB_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${LIBINT_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES 
    "${TiledArray_LINKER_FLAGS};${TiledArray_LIBRARIES}")

ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/molecule)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/basis)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/molecule)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/basis)

LIST(APPEND TileClusterChem_LIBRARIES molecule basis)

MESSAGE(STATUS "Using CXX Flags ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Using libraries ${TileClusterChem_LIBRARIES}")

ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/doc)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/tests)
enable_testing(true)

SET(TileClusterChem_Sources "main.cpp")

ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/examples)

ADD_EXECUTABLE(${PROJECT_NAME} ${TileClusterChem_Sources}) 
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
  ${TileClusterChem_LIBRARIES}
)
