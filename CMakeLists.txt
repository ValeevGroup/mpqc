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
# Compiler/library features
###############################################################################

include(CheckCXXSourceCompiles)
include(CheckIncludeFiles)
include(CheckFunctionExists)

check_include_files(sys/types.h          HAVE_SYS_TYPES_H)
check_include_files(sys/resource.h       HAVE_SYS_RESOURCE_H)
check_include_files(pwd.h                HAVE_PWD_H)
check_include_files(libunwind.h          HAVE_LIBUNWIND_H)
check_include_files(fenv.h               HAVE_FENV_H)
include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(cxxabi.h          HAVE_CXXABI_H)
if (HAVE_CXXABI_H)
  CHECK_CXX_SOURCE_COMPILES("
    #include <cxxabi.h>
    int main(int argc, char* argv[]) {
      char* result = abi::__cxa_demangle(\"\", 0, 0, 0);
      return 0;
    }" HAVE_CXA_DEMANGLE)
  if (HAVE_CXA_DEMANGLE)
    message("-- Looking for abi::__cxa_demangle - found")
  else()
    message("-- Looking for abi::__cxa_demangle - not found")
  endif()
endif(HAVE_CXXABI_H)
if (HAVE_LIBUNWIND_H)
  CHECK_FUNCTION_EXISTS (unw_init_local          HAVE_LIBUNWIND)
endif(HAVE_LIBUNWIND_H)
if (HAVE_EXECINFO_H)
  CHECK_FUNCTION_EXISTS (backtrace               HAVE_BACKTRACE_FN)
  CHECK_FUNCTION_EXISTS (backtrace_symbols       HAVE_BACKTRACE_SYMBOLS_FN)
  if (HAVE_BACKTRACE_FN AND HAVE_BACKTRACE_SYMBOLS_FN)
    set (HAVE_BACKTRACE)
  endif()
else (HAVE_EXECINFO_H)
  #find_package(ExecInfo)
endif (HAVE_EXECINFO_H)

CHECK_FUNCTION_EXISTS (gethostname       HAVE_GETHOSTNAME)
CHECK_FUNCTION_EXISTS (getpwuid          HAVE_GETPWUID)
CHECK_FUNCTION_EXISTS (geteuid           HAVE_GETEUID)
CHECK_FUNCTION_EXISTS (getrusage         HAVE_GETRUSAGE)
CHECK_FUNCTION_EXISTS (posix_spawn       HAVE_POSIX_SPAWN)
CHECK_FUNCTION_EXISTS (system            HAVE_SYSTEM)
CHECK_FUNCTION_EXISTS (signal            HAVE_SIGNAL)
CHECK_FUNCTION_EXISTS (feenableexcept    HAVE_FEENABLEEXCEPT)
CHECK_FUNCTION_EXISTS (fedisableexcept   HAVE_FEDISABLEEXCEPT)
CHECK_FUNCTION_EXISTS (setrlimit         HAVE_SETRLIMIT)

###############################################################################
# Dependencies
###############################################################################

# search for TA first, other packages can be found with it
FIND_PACKAGE(TiledArray CONFIG REQUIRED COMPONENTS tiledarray)
include(ValidateTiledArray)
validate_tiledarray()
set(TA_LIBRARIES tiledarray)
message(STATUS "Found TiledArray CONFIG at ${TiledArray_CONFIG}")

FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(Boost REQUIRED COMPONENTS serialization)
FIND_PACKAGE(Eigen3 REQUIRED)
FIND_PACKAGE(LIBINT2 REQUIRED)

INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LIBINT2_INCLUDE_DIRS})
add_definitions(${LIBINT2_EXTRA_DEFINITIONS})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

find_package(PythonInterp)

###############################################################################
# Declare LibMPQC
###############################################################################

INCLUDE_DIRECTORIES(${MPQC_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${MPQC_BINARY_DIR}/src)

ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/src/mpqc)

MESSAGE(STATUS "Using CXX Flags ${CMAKE_CXX_FLAGS}")

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
