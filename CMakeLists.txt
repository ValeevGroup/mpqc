cmake_minimum_required(VERSION 2.8.11)
project(TileClusterChem)

set(TileClusterChem_VERSION_MAJOR 0)
set(TileClusterChem_VERSION_MINOR 1)
SET(CMAKE_COLOR_MAKEFILE ON)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/cmake/modules/")

# search for TA first, other packages can be found with it
FIND_PACKAGE(TiledArray CONFIG HINTS ${TiledArray_CONFIG_DIR} REQUIRED)
INCLUDE_DIRECTORIES(${TiledArray_INCLUDE_DIRS})
LIST(APPEND TileClusterChem_LIBRARIES 
    "${TiledArray_LINKER_FLAGS};${TiledArray_LIBRARIES}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TiledArray_COMPILE_FLAGS}")

FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(Boost REQUIRED)
FIND_PACKAGE(LIBINT REQUIRED)
FIND_PACKAGE(BTAS REQUIRED)
FIND_PACKAGE(RapidJSON REQUIRED)

INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LIBINT_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${BTAS_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${RAPIDJSON_INCLUDE_DIRS})

ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CEREAL_INCLUDE_DIR})

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

LIST(APPEND TileClusterChem_LIBRARIES ${TBB_LIBRARIES})
LIST(APPEND TileClusterChem_LIBRARIES ${LIBINT_LIBRARIES})

ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/molecule)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/expression)
#ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/f12)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/integrals)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/basis)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/molecule)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/basis)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/expression)

LIST(APPEND TileClusterChem_LIBRARIES molecule basis integrals expression)

MESSAGE(STATUS "Using CXX Flags ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Using libraries ${TileClusterChem_LIBRARIES}")

ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/doc)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/tests)
enable_testing(true)

SET(TileClusterChem_Sources "main.cpp")

ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/examples)

ADD_EXECUTABLE(${PROJECT_NAME} ${TileClusterChem_Sources}) 
ADD_DEPENDENCIES(${PROJECT_NAME} cereal)
TARGET_LINK_LIBRARIES(${PROJECT_NAME}
  ${TileClusterChem_LIBRARIES}
)

# let clion know the source file
if(USING_CLION_AS_IDE)
#    include_directories(cc)
#    include_directories(integrals)
#    include_directories(mp2)
#    file(GLOB_RECURSE CLION_SRC "*.h" "*.cpp" "*.c" "*.cc" "*.hpp" "*.txt" ".in")
#    add_library(clion_get_sources EXCLUDE_FROM_ALL ${CLION_SRC})
#    set_target_properties(clion_get_sources PROPERTIES LINKER_LANGUAGE CXX)
    FILE (GLOB_RECURSE clion_all_headers ${CMAKE_SOURCE_DIR}/*.h)
    FILE (GLOB_RECURSE clion_all_source ${CMAKE_SOURCE_DIR}/*.cpp)
#    message(WARNING "${clion_all_headers}")
    ADD_CUSTOM_TARGET(all_clion
        SOURCES ${clion_all_headers} ${clion_all_source}
    )
    ADD_EXECUTABLE(clion_test "main.cpp" ${clion_all_headers})
    ADD_DEPENDENCIES(clion_test cereal)
    TARGET_LINK_LIBRARIES(clion_test ${TileClusterChem_LIBRARIES})
    #add_executable(dummy_name EXCLUDE_FROM_ALL ${clion_all_headers})
endif(USING_CLION_AS_IDE)
