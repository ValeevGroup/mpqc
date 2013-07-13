#; -*-CMake-*-

cmake_minimum_required (VERSION 2.8.8)
project (MPQC)

option(MPQC_EXPERT "MPQC Expert mode: disables automatically downloading or building dependencies" OFF)
option(MPQC_NEW_FEATURES "Use new MPQC features (requires C++11)" OFF)
option(MPQC_OPENMP "Use OpenMP" OFF)

enable_language (CXX)
if (NOT CMAKE_CXX_COMPILER)
  message(FATAL_ERROR "C++ compiler not found")
endif()

# FORTRAN setup
enable_language (Fortran)
if (NOT CMAKE_Fortran_COMPILER)
   message(FATAL_ERROR "Fortran compiler not found")
endif()
if(INTEGER8) 
  set(F77_INTEGER_WIDTH 8)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DINTEGER77='INTEGER*8' ")
else()
  set(F77_INTEGER_WIDTH 4)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DINTEGER77='INTEGER*4' ")
endif()

# extra cmake files are shipped with MPQC
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")

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

if (MPQC_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

include(TestCXXAcceptsFlag)
include(CheckCXXSourceCompiles)

if (MPQC_NEW_FEATURES)
  if (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    message("* C++11: assuming Clang doesn't need -std=c++0x flag.")
  else()
    # C++11 flag is required
    CHECK_CXX_ACCEPTS_FLAG("-std=c++0x" CXX_FLAG_CXX0X)
    CHECK_CXX_ACCEPTS_FLAG("-std=c++11" CXX_FLAG_CXX11)
    set(HAVE_CXX11 FALSE)
    if (CXX_FLAG_CXX0X)
      set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
      message(STATUS "C++11 flag: -std=c++0x")
      set(HAVE_CXX11 TRUE)
    elseif(CXX_FLAG_CXX11)
      set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
      message(STATUS "C++11 flag: -std=c++11")
      set(HAVE_CXX11 TRUE)
    endif()
  endif()
  add_definitions(-DMPQC_NEW_FEATURES)  
endif()

# C++11 features
CHECK_CXX_SOURCE_COMPILES("
  int main() {
  auto c = 1;
  }"
  HAVE_CXX_AUTO)  
if (MPQC_NEW_FEATURES AND NOT HAVE_CXX_AUTO)
  message(FATAL_ERROR "MPQC_NEW_FEATURES requires support for C++11 auto")
endif()

# std::array
CHECK_CXX_SOURCE_COMPILES("
    #include <array>
    int main(void)
    {
      std::array<int,2> array;
    }"
    HAVE_STD_ARRAY)

CHECK_CXX_SOURCE_COMPILES("
  int main() {
  void *restrict ptr;
  }" HAVE_CXX_RESTRICT)

CHECK_CXX_SOURCE_COMPILES("
  int main() {
  void *__restrict__ ptr;
  }" HAVE_CXX___RESTRICT__)



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
set(MPQC_MAJOR_VERSION 3)
set(MPQC_MINOR_VERSION 0)
set(MPQC_MICRO_VERSION 0)
set(MPQC_BUILDID alpha)

set(MPQC_VERSION "${MPQC_MAJOR_VERSION}.${MPQC_MINOR_VERSION}.${MPQC_MICRO_VERSION}")
set(MPQC_VERSION "${MPQC_VERSION}-${MPQC_BUILDID}")

include(CheckTypeSize)
check_type_size(int64_t C_TYPE_INT64_T)
check_type_size(int32_t C_TYPE_INT32_T)
check_type_size(long C_TYPE_LONG)
check_type_size(int C_TYPE_INT)

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

set(SRC_MPQC_DATA_PATH ${PROJECT_SOURCE_DIR}/lib)
set(MPQCDATAPATH ${CMAKE_INSTALL_FULL_DATAROOTDIR}/mpqc/${MPQC_VERSION})
set(INSTALLED_MPQC_DATA_PATH ${CMAKE_INSTALL_FULL_DATAROOTDIR}/mpqc/${MPQC_VERSION})

INSTALL(FILES ./lib/atominfo.kv DESTINATION ${MPQCDATAPATH})
INSTALL(DIRECTORY ./lib/basis DESTINATION ${MPQCDATAPATH})

include_directories("${PROJECT_BINARY_DIR}/src/lib")
include_directories("${PROJECT_SOURCE_DIR}/src/lib")
include_directories("${PROJECT_SOURCE_DIR}/include") # flex/lexer

# external dependencies
include(external/External)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/src/lib/mpqc_config.h.in
  ${PROJECT_BINARY_DIR}/src/lib/mpqc_config.h
)
add_definitions(-DHAVE_CONFIG_H)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/bin/mpqcrun.in
  ${PROJECT_BINARY_DIR}/bin/mpqcrun
  @ONLY
)

# libraries
set(LINK_LIBRARIES)
list(APPEND LINK_LIBRARIES ${LIBRYSQ_LIBRARY})
list(APPEND LINK_LIBRARIES ${LIBINT2_LIBRARY})
list(APPEND LINK_LIBRARIES ${MADNESS_LIBRARIES})
list(APPEND LINK_LIBRARIES ${Boost_LIBRARIES})
list(APPEND LINK_LIBRARIES ${LAPACK_LIBRARIES})
list(APPEND LINK_LIBRARIES ${BLAS_LIBRARIES})
list(APPEND LINK_LIBRARIES ${PYTHON_LIBRARIES})
list(APPEND LINK_LIBRARIES ${HDF5_LIBRARIES})
list(APPEND LINK_LIBRARIES ${MPI_LIBRARIES})
list(APPEND LINK_LIBRARIES ${PAPI_LIBRARIES})
list(APPEND LINK_LIBRARIES ${LIBRARIES})
list(APPEND LINK_LIBRARIES ${CMAKE_DL_LIBS})

# R-12
if (HAVE_LIBINT2)
  set(MPQC_R12 TRUE)
  add_definitions(-DMPQC_R12)
else()
  message("** R12 requires Libint2.  R12 methods are disabled")
endif()

# sources
add_subdirectory(src)
add_subdirectory(doc)

# checking/testing
add_subdirectory(test)


# Export include dirs, library list

# get_target_property(LINK_LIBRARIES mpqc LINK_LIBRARIES)
# message(STATUS "LIBRARIES = ${LIBRARIES}")

set(MPQC_CONFIG_LIBRARIES "")

# Fortran libraries
foreach (_lib ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  set(_libpath _libpath-NOTFOUND)
  find_library(_libpath ${_lib} PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  if (_libpath)
    list(APPEND LINK_LIBRARIES ${_libpath})
  else()
    list(APPEND LINK_LIBRARIES ${_lib})
  endif()
endforeach()

# transform library list into compiler args
foreach (_lib ${LINK_LIBRARIES})
  set(_libpath _libpath-NOTFOUND) # DON'T REMOVE THIS LINE
  find_library(_libpath ${_lib})
  #message(STATUS "_lib = ${_lib}, _libpath = ${_libpath}")
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
  set(MPQC_CONFIG_LIBRARIES "${MPQC_CONFIG_LIBRARIES} ${_library}")
endforeach()

# message(STATUS "MPQC_CONFIG_LIBRARIES = ${MPQC_CONFIG_LIBRARIES}")

# include dirs
get_directory_property(MPQC_INCLUDE_DIRECTORIES DIRECTORY src/bin/mpqc INCLUDE_DIRECTORIES)
# remove internal source dirs
list(REMOVE_ITEM
  MPQC_INCLUDE_DIRECTORIES
  ${PROJECT_BINARY_DIR}/src/lib
  ${PROJECT_SOURCE_DIR}/src/lib
  ${PROJECT_SOURCE_DIR}/include
) 

set(MPQC_CPP_FLAGS "${CMAKE_CPP_FLAGS}")
foreach (_dir ${MPQC_INCLUDE_DIRECTORIES})
  set(MPQC_CPP_FLAGS "${MPQC_CPP_FLAGS} -I${_dir}")
endforeach()
# MPQC includes
set(MPQC_CPP_FLAGS "-I${CMAKE_INSTALL_PREFIX}/include/ ${MPQC_CPP_FLAGS}")
set(MPQC_CPP_FLAGS "-I${CMAKE_INSTALL_PREFIX}/include/mpqc ${MPQC_CPP_FLAGS}")
#message(STATUS "MPQC_CPP_FLAGS = ${MPQC_CPP_FLAGS}")

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/bin/mpqc-config.in
  ${PROJECT_BINARY_DIR}/bin/mpqc-config
  @ONLY
)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/mpqc-config DESTINATION bin)
