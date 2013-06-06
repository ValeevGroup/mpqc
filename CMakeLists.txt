#; -*-CMake-*-

enable_language (CXX)
enable_language (Fortran)

cmake_minimum_required (VERSION 2.8.8)
project (MPQC)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-g ${CMAKE_CXX_FLAGS_RELEASE}")

set(CMAKE_SKIP_RPATH FALSE)

set(BUILD_TESTING FALSE CACHE BOOLEAN "BUILD_TESTING")
set(BUILD_TESTING_STATIC FALSE CACHE BOOLEAN "BUILD_TESTING_STATIC")
set(BUILD_TESTING_SHARED FALSE CACHE BOOLEAN "BUILD_TESTING_SHARED")

# string(TOUPPER ${CMAKE_SYSTEM_NAME} TARGET_)
# if(ARCH EQUAL 64)
#   set(TARGET_ "${TARGET_}64")
# endif()

# set(TARGET "${TARGET_}" CACHE STRING "TARGET")
# add_definitions(-D${TARGET})
# if("${TARGET}" MATCHES "64")
#   string(REGEX REPLACE "64" "" TARGET_ "${TARGET}")
#   add_definitions(-D${TARGET_})
# endif()
# MESSAGE(STATUS "TARGET: ${TARGET}")

# if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#  set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
# endif()

option(OPENMP "Enable OpenMP" ON)
if(OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()


include(CheckCXXSourceCompiles)
include(TestCXXAcceptsFlag)

CHECK_CXX_SOURCE_COMPILES("
  int main() {
  void *restrict ptr;
  }" HAVE_CXX_RESTRICT)

CHECK_CXX_SOURCE_COMPILES("
  int main() {
  void *__restrict__ ptr;
  }" HAVE_CXX___RESTRICT__)


option(MPQC_NEW_FEATURES "Use new MPQC features (requires C++11)" OFF)
if (MPQC_NEW_FEATURES)

  add_definitions(-DMPQC_NEW_FEATURES)

  # C++11 flag
  CHECK_CXX_ACCEPTS_FLAG("-std=c++0x" CXX_FLAG_CXX0X)
  CHECK_CXX_ACCEPTS_FLAG("-std=c++11" CXX_FLAG_CXX11)

  if (CXX_FLAG_CXX0X)
    set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
    MESSAGE(STATUS "C++11 flag: -std=c++0x")
  elseif(CXX_FLAG_CXX11)
    set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    MESSAGE(STATUS "C++11 flag: -std=c++11")
  endif()

  # C++11 features
  CHECK_CXX_SOURCE_COMPILES("
    int main() {
    auto c = 1;
    }" HAVE_CXX_AUTO)
  if (NOT HAVE_CXX_AUTO)
    message(FATAL_ERROR "No support for C++ auto")
  endif()

  # rvalue references
  CHECK_CXX_SOURCE_COMPILES("
    #include <utility>
    int foo(int &&i) { return i; }
    int main(void)
    {
      return foo(std::forward<int>(0));
    }" HAVE_CXX_RVALUE
    )
  if (NOT HAVE_CXX_RVALUE)
    message(FATAL_ERROR "No support for C++ rvalue/std::forward")
  endif()

endif(MPQC_NEW_FEATURES)


option(USE_STATIC_LIBS "use static libraries" OFF)
if(USE_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
ENDIF()

# find_package(PythonLibs)
# if(PYTHONLIBS_FOUND)
#   set(HAVE_PYTHON TRUE)
#   include_directories(${PYTHON_INCLUDE_DIR})
# endif()


find_package(Threads)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(LIBRARIES "${LIBRARIES};${CMAKE_THREAD_LIBS_INIT}")
endif()


set(TARGET_ARCH "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
set(SC_MAJOR_VERSION 3)
set(SC_MINOR_VERSION 0)
set(SC_MICRO_VERSION 0)
set(SC_BUILDID alpha)

set(SC_VERSION "${SC_MAJOR_VERSION}.${SC_MINOR_VERSION}.${SC_MICRO_VERSION}")
set(SC_VERSION "${SC_VERSION}-${SC_BUILDID}")

check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(dlfcn.h HAVE_DLFCN_H)
check_include_files(sys/time.h HAVE_SYS_TIME_H)
check_include_files(sys/times.h HAVE_SYS_TIMES_H)
check_include_files(sys/resource.h HAVE_SYS_RESOURCE_H)
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(pwd.h HAVE_PWD_H)
check_include_files(time.h HAVE_TIME_H)

include(CheckFunctionExists)

CHECK_FUNCTION_EXISTS (gethostname       HAVE_GETHOSTNAME)
CHECK_FUNCTION_EXISTS (getpwuid          HAVE_GETPWUID)
CHECK_FUNCTION_EXISTS (geteuid           HAVE_GETEUID)
CHECK_FUNCTION_EXISTS (getrusage         HAVE_GETRUSAGE)
CHECK_FUNCTION_EXISTS (lstat             HAVE_LSTAT)
CHECK_FUNCTION_EXISTS (time              HAVE_TIME)
CHECK_FUNCTION_EXISTS (ctime             HAVE_CTIME)

set(BLAS_F77_INTEGER_WIDTH 4)
if(INTEGER8) 
  set(BLAS_F77_INTEGER_WIDTH 8)
endif()  

include(FortranCInterface)

set(FORTRAN_FUNCTIONS
  PDSTEQR
  DCOPY
  DNRM2
  DSCAL
  DGEMM
  DGEMV
  DAXPY
  DDOT
  DSPMV
  DGESVD
  DSPSVX
  DSYEVD
  DSPTRF
  DPPTRF
  DSPTRI
  DPPTRI
  DLANSP
  DSPCON
  DPPCON
  DLAMCH
  DLACPY
  DSPTRS
  DPPTRS
  DSPRFS
  DPPRFS
  DSYGV
  )

set(fortran_headers 
lib/math/optimize/f77sym.h
lib/math/optimize/levmar/f77sym.h
lib/math/scmat/f77sym.h
)

foreach(header ${fortran_headers})
  FortranCInterface_HEADER(
    ${PROJECT_SOURCE_DIR}/src/${header}
    MACRO_NAMESPACE "FC_"
    SYMBOL_NAMESPACE "F77_"
    SYMBOLS ${FORTRAN_FUNCTIONS}
    )
endforeach()

include(GNUInstallDirs)           

set(SRC_SCLIBDIR ${PROJECT_SOURCE_DIR}/src/lib)
set(SCDATADIR ${CMAKE_INSTALL_FULL_DATAROOTDIR}/mpqc/${SC_VERSION})
set(INSTALLED_SCLIBDIR ${CMAKE_INSTALL_FULL_LIBDIR})

INSTALL(FILES ./lib/atominfo.kv DESTINATION ${SCDATADIR})
INSTALL(DIRECTORY ./lib/basis DESTINATION ${SCDATADIR})


include_directories("${PROJECT_SOURCE_DIR}/include") # flex/lexer
include_directories("${PROJECT_SOURCE_DIR}/src/lib")
include_directories("${PROJECT_BINARY_DIR}/src/lib")


find_package(MPI)
if (NOT MPI_FOUND)
  set(MPI_LIBRARIES "")
  MESSAGE(ERROR "MPI not found")
else()
  set(HAVE_MPI TRUE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINKER_FLAGS}")
  set(LIBRARIES "${LIBRARIES};${MPI_LIBRARIES}")
  include_directories(${MPI_INCLUDE_PATH})
  message(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
  message(STATUS "MPI_LINKER_FLAGS: ${MPI_LINKER_FLAGS}")
  #add_definitions(-DOMPI_SKIP_MPICXX)
  option(ALWAYS_USE_MPI "Always use MPI" TRUE)
endif()

# external dependencies
include(external/CMakeLists.txt)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/src/lib/scconfig.h.in
  ${PROJECT_BINARY_DIR}/src/lib/scconfig.h
)
add_definitions(-DHAVE_CONFIG_H)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/src/lib/scdirlist.h.in
  ${PROJECT_BINARY_DIR}/src/lib/scdirlist.h
)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/bin/mpqcrun.in
  ${PROJECT_BINARY_DIR}/bin/mpqcrun
  @ONLY
)

# # ?is the LocalMakefile necessary
# CONFIGURE_FILE(
#   ${PROJECT_SOURCE_DIR}/lib/LocalMakefile.in
#   ${PROJECT_BINARY_DIR}/lib/LocalMakefile
# )


add_subdirectory(src)

# checking/testing
add_subdirectory(test)
