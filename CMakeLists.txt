cmake_minimum_required(VERSION 3.0.0)

###############################################################################
# Announce ourselves
###############################################################################
project(MPQC)

set (MPQC_MAJOR_VERSION 4)
set (MPQC_MINOR_VERSION 0)
set (MPQC_PATCH_VERSION 0)
set(MPQC_BUILDID alpha)

set(MPQC_VERSION "${MPQC_MAJOR_VERSION}.${MPQC_MINOR_VERSION}.${MPQC_PATCH_VERSION}")
if (MPQC_BUILDID)
  set(MPQC_VERSION "${MPQC_VERSION}-${MPQC_BUILDID}")
endif(MPQC_BUILDID)

set(TARGET_ARCH "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

###############################################################################
# Load up extra CMake packages
###############################################################################
include(FeatureSummary)

###############################################################################
# Options
###############################################################################
# boolean
# non-boolean
if (NOT DEFINED MPQC_MEMORY_CHECK)
  set(MPQC_MEMORY_CHECK -1) # How to check memory use (see mpqc_config.h for details)
else()
  message("MPQC_MEMORY_CHECK = ${MPQC_MEMORY_CHECK}")
endif()
if (NOT DEFINED MPQC_ASSERT_MODE)
  set(MPQC_ASSERT_MODE -1) # How to implement assertions (see mpqc_config.h for details)
endif()

###############################################################################
# How to compile
###############################################################################

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

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

###############################################################################
# Dependencies
###############################################################################

# grand list of all MPQC libraries and dependencies goes here; may go away
set(MPQC_LIBRARIES "" CACHE INTERNAL "MPQC libs and dependencies")

# search for TA first, other packages can be found with it
FIND_PACKAGE(TiledArray CONFIG REQUIRED COMPONENTS tiledarray)
set(TA_LIBRARIES tiledarray)
set(MPQC_LIBRARIES "${TA_LIBRARIES}")
message(STATUS "Found TiledArray CONFIG at ${TiledArray_CONFIG}")

FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS serialization)
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(LIBINT2 REQUIRED)
FIND_PACKAGE(RapidJSON)

INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LIBINT2_INCLUDE_DIRS})
add_definitions(${LIBINT2_EXTRA_DEFINITIONS})
INCLUDE_DIRECTORIES(${RAPIDJSON_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

set(MPQC_LIBRARIES
    "${MPQC_LIBRARIES};${TBB_LIBRARIES};${LIBINT2_LIBRARIES};${Boost_LIBRARIES}")

###############################################################################
# Declare LibMPQC
###############################################################################

INCLUDE_DIRECTORIES(${MPQC_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${MPQC_BINARY_DIR}/src)

ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/src/mpqc)

LIST(APPEND MPQC_LIBRARIES MPQCmolecule MPQCbasis MPQCwfn MPQCproperties MPQCmath_linalg MPQCmath_ta MPQCmath_clr MPQCexpression MPQCf12 MPQCcc MPQCmbpt MPQCintegrals MPQCscf MPQCkeyval)

message (STATUS "MPQC_LIBRARIES=${MPQC_LIBRARIES}")

MESSAGE(STATUS "Using CXX Flags ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "Using libraries ${MPQC_LIBRARIES}")

###############################################################################
# Declare MPQC program
###############################################################################

# none yet
add_subdirectory(examples)
ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/src/bin)

###############################################################################
# Declare MPQC tests and doc
###############################################################################

ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/doc)
ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/tests)
enable_testing(true)

###############################################################################
# Wrap up
###############################################################################

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/src/mpqc/mpqc_config.h.in
  ${PROJECT_BINARY_DIR}/src/mpqc/mpqc_config.h
)

feature_summary(WHAT ALL
                DESCRIPTION "=== MPQC Package/Feature Info ===")
                
###############################################################################
# appendix: misc details
###############################################################################
SET(CMAKE_COLOR_MAKEFILE ON)

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
    #ADD_EXECUTABLE(clion_test "main.cpp" ${clion_all_headers})
    #ADD_DEPENDENCIES(clion_test cereal)
    #TARGET_LINK_LIBRARIES(clion_test ${MPQC_LIBRARIES})
    #add_executable(dummy_name EXCLUDE_FROM_ALL ${clion_all_headers})
endif(USING_CLION_AS_IDE)
