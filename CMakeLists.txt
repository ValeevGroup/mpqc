#; -*-CMake-*-

enable_language (CXX)

enable_language (Fortran)
if (NOT CMAKE_Fortran_COMPILER)
   message(FATAL_ERROR "Fortran compiler not found")
endif()

cmake_minimum_required (VERSION 2.8.8)
project (MPQC)

# extra cmake files are shipped with MPQC
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-g ${CMAKE_CXX_FLAGS_RELEASE}")

set(CMAKE_SKIP_RPATH FALSE)

set(BUILD_TESTING FALSE CACHE BOOLEAN "BUILD_TESTING")
set(BUILD_TESTING_STATIC FALSE CACHE BOOLEAN "BUILD_TESTING_STATIC")
set(BUILD_TESTING_SHARED FALSE CACHE BOOLEAN "BUILD_TESTING_SHARED")

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

endif(MPQC_NEW_FEATURES)

# C++11 features
CHECK_CXX_SOURCE_COMPILES("
    int main() {
    auto c = 1;
    }"
    HAVE_CXX_AUTO)  

# rvalue references
CHECK_CXX_SOURCE_COMPILES("
    #include <utility>
    int foo(int &&i) { return i; }
    int main(void)
    {
      return foo(std::forward<int>(0));
    }"
    HAVE_CXX_RVALUE)

# std::array
CHECK_CXX_SOURCE_COMPILES("
    #include <array>
    int main(void)
    {
      std::array<int,2> array;
    }"
    HAVE_STD_ARRAY)

if (MPQC_NEW_FEATURES)
  if (NOT HAVE_CXX_AUTO)
    message(FATAL_ERROR "No support for C++ auto")
  endif()
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
if (Threads_FOUND)
  #message(STATUS "Threads: ${CMAKE_THREAD_LIBS_INIT}")
  set(LIBRARIES "${LIBRARIES};${CMAKE_THREAD_LIBS_INIT}")
  check_include_files(pthread.h HAVE_PTHREAD)
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
FortranCInterface_VERIFY(CXX)

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
  DGEQRF
  DGETRF
  DGETRS
  DORGQR
  DPOTF2
  DTRTRS
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

include_directories("${PROJECT_BINARY_DIR}/src/lib")
include_directories("${PROJECT_SOURCE_DIR}/src/lib")
include_directories("${PROJECT_SOURCE_DIR}/include") # flex/lexer

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
include(external/External)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/src/lib/scconfig.h.in
  ${PROJECT_BINARY_DIR}/src/lib/scconfig.h
)
add_definitions(-DHAVE_CONFIG_H)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/bin/mpqcrun.in
  ${PROJECT_BINARY_DIR}/bin/mpqcrun
  @ONLY
)

# libraries
set(LINK_LIBRARIES "")
if (HAVE_LIBINT2)
  list(APPEND LINK_LIBRARIES ${LIBINT2_LIBRARY})
endif()
list(APPEND LINK_LIBRARIES ${MADNESS_LIBRARIES})
list(APPEND LINK_LIBRARIES ${Boost_LIBRARIES})
list(APPEND LINK_LIBRARIES ${LAPACK_LIBRARIES})
list(APPEND LINK_LIBRARIES ${BLAS_LIBRARIES})
list(APPEND LINK_LIBRARIES ${PYTHON_LIBRARIES})
list(APPEND LINK_LIBRARIES ${LIBRARIES})
list(APPEND LINK_LIBRARIES ${CMAKE_DL_LIBS})

# sources
add_subdirectory(src)
add_subdirectory(doc)

# checking/testing
add_subdirectory(test)


# Export include dirs, library list

# get_target_property(LINK_LIBRARIES mpqc LINK_LIBRARIES)
# message(STATUS "LINK_LIBRARIES = ${LINK_LIBRARIES}")

set(MPQC_LINK_LIBRARIES "")

# Fortran libraries
foreach (_lib ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  set(_libpath _path-NOTFOUND)
  find_library(_libpath ${_lib} PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  list(APPEND LINK_LIBRARIES ${_libpath})
endforeach()

# transform library list into compiler args
foreach (_lib ${LINK_LIBRARIES})
  set(_libpath _path-NOTFOUND)
  find_library(_libpath ${_lib} PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  #message(STATUS "_lib = ${_lib}, _path = ${_path}")
  if (_libpath)
    set(_library "${_libpath}")
  else()
    set(_library "${_lib}")
  endif()
  # Apple framework libs
  get_filename_component(_ext ${_library} EXT)
  if ("${_ext}" STREQUAL ".framework")
    get_filename_component(_name ${_library} NAME_WE)
    get_filename_component(_path ${_library} PATH)
    #message(STATUS "${_library} = ${_name} ${_path} ${_ext}")
    set(_library "-F ${_path} -framework ${_name}")
  endif()
  #message(STATUS "${_lib} = ${_library}")
  set(MPQC_LINK_LIBRARIES "${MPQC_LINK_LIBRARIES} ${_library}")
endforeach()

# message(STATUS "mpqc_link_libraries = ${MPQC_LINK_LIBRARIES}")

# include dirs sans the source dirs
get_directory_property(MPQC_INCLUDE_DIRECTORIES DIRECTORY src/bin/mpqc INCLUDE_DIRECTORIES)
list(REMOVE_ITEM
  MPQC_INCLUDE_DIRECTORIES
  ${PROJECT_BINARY_DIR}/src/lib
  ${PROJECT_SOURCE_DIR}/src/lib
  ${PROJECT_SOURCE_DIR}/include
) 

set(MPQC_CPP_FLAGS "${CMAKE_CPP_FLAGS}")
foreach (_dir ${MPQC_INCLUDE_DIRECTORIES} ${CMAKE_INSTALL_PREFIX}/include/mpqc)
  set(MPQC_CPP_FLAGS "${MPQC_CPP_FLAGS} -I${_dir}")
endforeach()
#message(STATUS "MPQC_CPP_FLAGS = ${MPQC_CPP_FLAGS}")

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/bin/sc-config.in
  ${PROJECT_BINARY_DIR}/bin/sc-config
  @ONLY
)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/sc-config DESTINATION bin)
