cmake_minimum_required(VERSION 3.0.0)
project(TileClusterChem)

set(TileClusterChem_VERSION_MAJOR 0)
set(TileClusterChem_VERSION_MINOR 1)
SET(CMAKE_COLOR_MAKEFILE ON)

if(BUILD_SHARED_LIBS)
  set(BLA_STATIC FALSE)
  set(CMAKE_MACOSX_RPATH TRUE)
else()
  set(BLA_STATIC TRUE)
  set(CMAKE_MACOSX_RPATH FALSE)
endif()

set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG "-Wl,-rpath,")
  # look for frameworks and appbundles last
  set(CMAKE_FIND_FRAMEWORK LAST)
  set(CMAKE_FIND_APPBUNDLE LAST)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/cmake/modules/")

# search for TA first, other packages can be found with it
FIND_PACKAGE(TiledArray CONFIG REQUIRED COMPONENTS tiledarray)
set(TA_LIBRARIES tiledarray)
LIST(APPEND TileClusterChem_LIBRARIES ${TA_LIBRARIES})

FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(Boost REQUIRED)
FIND_PACKAGE(LIBINT REQUIRED)
FIND_PACKAGE(RapidJSON REQUIRED)

INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LIBINT_INCLUDE_DIRS})
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
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/f12)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/integrals)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/basis)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/tensor)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/utility)
ADD_SUBDIRECTORY(${TileClusterChem_SOURCE_DIR}/scf)

LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/tensor)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/molecule)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/utility)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/basis)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/expression)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/f12)
LINK_DIRECTORIES(${TileClusterChem_BINARY_DIR}/scf)

LIST(APPEND TileClusterChem_LIBRARIES molecule basis tensor utility f12 integrals expression scf)

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
