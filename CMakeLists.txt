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

set(TARGET_ARCH "${CMAKE_SYSTEM}-${CMAKE_SYSTEM_PROCESSOR}")

###############################################################################
# Load up extra CMake packages
###############################################################################
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules/")
include(FeatureSummary)
include(GNUInstallDirs)
include(AddMPQCLibrary)
include(CMakePackageConfigHelpers)
include(AppendFlags)
#include(cotire)

# Set install paths ============================================================

set(MPQC_INSTALL_BINDIR "bin" 
    CACHE PATH "MPQC binary install directory")
set(MPQC_INSTALL_INCLUDEDIR "include" 
    CACHE PATH "MPQC INCLUDE install directory")
set(MPQC_INSTALL_LIBDIR "lib"
    CACHE PATH "MPQC LIB install directory")
set(MPQC_INSTALL_DATADIR "share/mpqc/${MPQC_VERSION}/data"
    CACHE PATH "MPQC DATA install directory")
# TODO this matches Automake value, should be adjusted to the convention: share/mpqc/${MPQC_VERSION}/doc
set(MPQC_INSTALL_DOCDIR "share/doc/mpqc"
    CACHE PATH "MPQC DOC install directory")
set(MPQC_INSTALL_CMAKEDIR "lib/cmake/mpqc"
    CACHE PATH "MPQC CMAKE install directory")

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

set(TA_POLICY sparse CACHE STRING "Set to \"dense\" (to default to TA::DensePolicy) or \"sparse\" (to default to TA::SparsePolicy)")
# convert string value of TA_POLICY to numerical values as TA_DEFAULT_POLICY
set(TA_DEFAULT_POLICY 1)
if(TA_POLICY STREQUAL sparse)
  set(TA_DEFAULT_POLICY 1)
  MESSAGE("Set to TA::SparsePolicy")
elseif(TA_POLICY STREQUAL dense)
  set(TA_DEFAULT_POLICY 0)
  MESSAGE("Set to TA::DensePolicy")
endif()
##########################
# Standard build variables
##########################

# Get standard build variables from the environment if they have not already been set
if(NOT CMAKE_C_FLAGS OR NOT DEFINED CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS "$ENV{CPPFLAGS}")
  append_flags(CMAKE_C_FLAGS "$ENV{CFLAGS}")
endif()
if(NOT CMAKE_CXX_FLAGS OR NOT DEFINED CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "$ENV{CPPFLAGS}")
  append_flags(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS}")
endif()
if(NOT CMAKE_EXE_LINKER_FLAGS OR NOT DEFINED CMAKE_EXE_LINKER_FLAGS)
  set(CMAKE_EXE_LINKER_FLAGS "$ENV{LDFLAGS}")
endif()

enable_language (CXX)
if (NOT CMAKE_CXX_COMPILER)
  message(FATAL_ERROR "C++ compiler not found")
endif()

###############################################################################
# How to compile
###############################################################################

### RPATH config to ensure full path is used
### see https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH
set(CMAKE_MACOSX_RPATH FALSE)
# use RPATH if possible
set(CMAKE_SKIP_RPATH FALSE)
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# MacOS specific
if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set(CMAKE_SHARED_LIBRARY_RUNTIME_C_FLAG "-Wl,-rpath,")
  # look for frameworks and appbundles last
  set(CMAKE_FIND_FRAMEWORK LAST)
  set(CMAKE_FIND_APPBUNDLE LAST)
endif()

# may need this to link in MPI library via MADWorld dependency alone
# e.g. when MADNESS is configured with non-throwing asserts there is no direct mention of MPI
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  set (CMAKE_EXE_LINKER_FLAGS "-Wl,--no-as-needed")
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# by default, use same type of BLAS/LAPACK libs (shared/static) as we are using for MPQC
if(BUILD_SHARED_LIBS)
  set(BLA_STATIC FALSE)
else()
  set(BLA_STATIC TRUE)
endif()

# look for ccache
find_program(CCACHE ccache)
if(CCACHE)
    message (STATUS "Found ccache: ${CCACHE}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE})
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ${CCACHE})
endif(CCACHE)

###############################################################################
# Compiler/library features
###############################################################################

include(CheckCXXSourceCompiles)
include(CheckIncludeFiles)
include(CheckFunctionExists)

# Check compiler feature support ===============================================

# need for C++14 will be inherited from madness
# but insist on strict standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

check_include_files(sys/types.h          HAVE_SYS_TYPES_H)
check_include_files(sys/resource.h       HAVE_SYS_RESOURCE_H)
check_include_files(pwd.h                HAVE_PWD_H)
check_include_files(libunwind.h          HAVE_LIBUNWIND_H)
check_include_files(execinfo.h           HAVE_EXECINFO_H)
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
    set (HAVE_BACKTRACE 1)
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

##### TiledArray #####
# search for TA first, other packages can be found with it
if(TiledArray_INSTALL_DIR)
  set(TiledArray_DIR ${TiledArray_INSTALL_DIR}/lib/cmake/tiledarray)
endif()
FIND_PACKAGE(TiledArray CONFIG REQUIRED COMPONENTS tiledarray)
include(ValidateTiledArray)
validate_tiledarray()
message(STATUS "Found TiledArray CONFIG at ${TiledArray_CONFIG}")

##### Boost #####
FIND_PACKAGE(Boost REQUIRED COMPONENTS filesystem locale serialization system)

##### TBB #####
# TiledArray and MADNESS set search path for TBB to TBB_ROOT_DIR, KitWare's FindTBB expects TBB_ROOT
if(TBB_ROOT_DIR)
  set(TBB_ROOT ${TBB_ROOT_DIR})
endif(TBB_ROOT_DIR)
FIND_PACKAGE(TBB REQUIRED)
# set up an interface library for TBB a la https://github.com/justusc/FindTBB
include(ImportTBB)
import_tbb()

##### Eigen3 #####
include(external/eigen.cmake)

##### BTAS #####
include(external/btas.cmake)

##### Libint2 #####
FIND_PACKAGE(LIBINT2 REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${LIBINT2_INCLUDE_DIRS})
add_definitions(${LIBINT2_EXTRA_DEFINITIONS})

find_package(PythonInterp)

# pass on some info to pkg-config
list(APPEND MPQC_CONFIG_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${LIBINT2_INCLUDE_DIRS})
list(APPEND MPQC_CONFIG_LIBRARIES ${Boost_LIBRARIES})
list(APPEND MPQC_CONFIG_DEFINITIONS ${LIBINT2_EXTRA_DEFINITIONS})

###############################################################################
# Descend into the source tree and build component libraries and MPQC exec
###############################################################################

INCLUDE_DIRECTORIES(${MPQC_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${MPQC_BINARY_DIR}/src)
ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/src/mpqc)

MESSAGE(STATUS "Using CXX Flags ${CMAKE_CXX_FLAGS}")

###############################################################################
# Declare MPQC program
###############################################################################

ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/src/bin)

###############################################################################
# Declare MPQC tests, examples, and docs
###############################################################################

ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/doc)
ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/tests)
ADD_SUBDIRECTORY(${MPQC_SOURCE_DIR}/examples)

##########################
# pkg-config variables
##########################
foreach(_inc ${MPQC_CONFIG_INCLUDE_DIRS})
  append_flags(MPQC_PC_CFLAGS "-I${_inc}")
endforeach()
foreach(_def ${MPQC_CONFIG_DEFINITIONS})
  append_flags(MPQC_PC_CFLAGS "${_def}")
endforeach()
foreach(_lib ${MPQC_CONFIG_LIBRARIES})
  append_flags(MPQC_PC_LIBS "${_lib}")
endforeach()

##########################
# Get the git revision tag information
##########################

if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
  find_package(Git REQUIRED)
  execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse -q HEAD
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE MPQC_REVISION )
  string(REGEX MATCH "[0-9a-f]*"
      MPQC_REVISION "${MPQC_REVISION}")
else()
  set(MPQC_REVISION "${MPQC_REVISION}")
endif()

###############################################################################
# Wrap up
###############################################################################

# Force cache refresh for compile flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "C compile flags" FORCE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "C++ compile flags" FORCE)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/src/mpqc/mpqc_config.h.in
  ${PROJECT_BINARY_DIR}/src/mpqc/mpqc_config.h
)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/mpqc.pc.in
  ${PROJECT_BINARY_DIR}/mpqc.pc
)

# sample CMakeLists.txt
set(sample_dirs "mp2")
foreach(sample_dir ${sample_dirs})
  configure_file(
    ${PROJECT_SOURCE_DIR}/doc/devsamp/${sample_dir}/CMakeLists.txt.in
    ${PROJECT_BINARY_DIR}/doc/devsamp/${sample_dir}/CMakeLists.txt
    @ONLY
  )
endforeach()

# install config files
install(FILES ${PROJECT_BINARY_DIR}/mpqc.pc
    DESTINATION lib/pkgconfig)

# Create the version file
write_basic_package_version_file(mpqc-config-version.cmake
  VERSION ${MPQC_VERSION} COMPATIBILITY AnyNewerVersion)

# Create the targets file
export(EXPORT mpqc
  FILE "${PROJECT_BINARY_DIR}/mpqc-targets.cmake")

# Create the configure file
configure_package_config_file(cmake/mpqc-config.cmake.in
    "${PROJECT_BINARY_DIR}/mpqc-config.cmake"
  INSTALL_DESTINATION "${MPQC_INSTALL_CMAKEDIR}"
  PATH_VARS CMAKE_INSTALL_PREFIX MPQC_INSTALL_BINDIR 
            MPQC_INSTALL_INCLUDEDIR MPQC_INSTALL_LIBDIR
            MPQC_INSTALL_DOCDIR MPQC_INSTALL_CMAKEDIR)

# Install config, version, and target files
install(EXPORT mpqc
    FILE "mpqc-targets.cmake"
    DESTINATION "${MPQC_INSTALL_CMAKEDIR}" 
    COMPONENT mpqc-config)
install(FILES
    "${PROJECT_BINARY_DIR}/mpqc-config.cmake"
    "${PROJECT_BINARY_DIR}/mpqc-config-version.cmake"
    DESTINATION "${MPQC_INSTALL_CMAKEDIR}" 
    COMPONENT mpqc-config)
add_custom_target(install-config
     COMMAND ${CMAKE_COMMAND} -DCOMPONENT=mpqc-config -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
     COMMENT "Installing MPQC config components")

feature_summary(WHAT ALL
                DESCRIPTION "=== MPQC Package/Feature Info ===")
                
###############################################################################
# appendix: misc details
###############################################################################
SET(CMAKE_COLOR_MAKEFILE ON)
