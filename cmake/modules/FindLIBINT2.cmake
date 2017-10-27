
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

if (NOT ${EIGEN3_FOUND})
  find(Eigen3 REQUIRED)
endif()

if (LIBINT2_INCLUDE_DIRS)
  set (LIBINT2_FOUND TRUE)

else (LIBINT2_INCLUDE_DIRS)
  set (LIBINT2_FOUND FALSE)

  if(LIBINT2_INSTALL_DIR)
    set(_LIBINT2_INSTALL_DIR ${LIBINT2_INSTALL_DIR})
  endif(LIBINT2_INSTALL_DIR)

  set(LIBINT2_INC_SEARCH_DIR ${_LIBINT2_INSTALL_DIR}/include)

  find_path(LIBINT2_INCLUDE_DIR NAMES libint2.hpp
    HINTS
    ${LIBINT2_INC_SEARCH_DIR})

  mark_as_advanced(LIBINT2_INCLUDE_DIR)

  # force linking against shared Libint only if building shared MPQC libs
  # to avoid linking errors (best case) or multiple instances of libint global data (worst case ... looking at you MacOS)
  if (BUILD_SHARED_LIBS)
    message(WARNING "Since BUILD_SHARED_LIBS=TRUE, will only look for SHARED libint2 library (static library will be ignored)")
    set(_LIBINT2_LIB_NAMES "libint2.so" "libint2.dylib")
  else (BUILD_SHARED_LIBS)
    set(_LIBINT2_LIB_NAMES "int2")
  endif(BUILD_SHARED_LIBS)

  set(_LIBINT2_LIBRARY_DIR ${_LIBINT2_INSTALL_DIR}/lib)
  set(LIBINT2_LIBRARY "")

  set(current_lib NOTFOUND)
  find_library(current_lib NAMES ${_LIBINT2_LIB_NAMES} HINTS ${_LIBINT2_LIBRARY_DIR})
  if(NOT current_lib)
    message(FATAL_ERROR "MPQC could not find Libint2 lib: tried ${_LIBINT2_LIB_NAMES}")
  endif(NOT current_lib)
  list(APPEND LIBINT2_LIBRARY ${current_lib})

  set(LIBINT2_LIBRARY_DIR ${_LIBINT2_LIBRARY_DIR})
  mark_as_advanced(LIBINT2_LIBRARY)

  # CODA
  if(LIBINT2_INCLUDE_DIR AND LIBINT2_LIBRARY)
  
    # validate version, etc. by compiling tests
    cmake_push_check_state()
  
    list(APPEND CMAKE_REQUIRED_INCLUDES ${EIGEN3_INCLUDE_DIR} ${LIBINT2_INCLUDE_DIR})
    list(APPEND CMAKE_REQUIRED_LIBRARIES ${LIBINT2_LIBRARY})
    list(APPEND CMAKE_REQUIRED_DEFINITIONS ${CMAKE_CXX14_STANDARD_COMPILE_OPTION})
    
    # sanity check
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <libint2.hpp>
      int main(int argc, char** argv) {
        libint2::initialize();
        libint2::finalize();
        return 0;
      }
      "  LIBINT2_COMPILES)
    if (NOT LIBINT2_COMPILES)
      message(FATAL_ERROR "Could not compile LIBINT2 test program.\nSee CMakeFiles/CMakeError.log for details")
    endif()
    
    # make sure libint2 version is up to date
    CHECK_CXX_SOURCE_COMPILES(
      "
          #include <libint2.hpp>
          #include <libint2/engine.h>
          #include <libint2/util/any.h>
          #if !((LIBINT_MAJOR_VERSION>2) || (LIBINT_MAJOR_VERSION==2 && LIBINT_MINOR_VERSION>4) || (LIBINT_MAJOR_VERSION==2 && LIBINT_MINOR_VERSION==4 && LIBINT_MICRO_VERSION>=0))
          # error \"Libint2 library is too old\"
          #endif
          int main(int argc, char** argv) {
          
            // need libint2::Engine
            libint2::Engine e0;
            
            // check for standard-conforming libint2::any
            libint2::any x = 0;
            bool has_value_exists = x.has_value();
            
            return 0;
          }
      "  LIBINT2_IS_UP_TO_DATE)    
    if (NOT LIBINT2_IS_UP_TO_DATE)
      message(FATAL_ERROR "Libint2 library is too old: 2.4.0 is required")
    endif()

    # make sure libint2 is properly configured
    CHECK_CXX_SOURCE_COMPILES(
      "
          #include <libint2.hpp>
          #include <libint2/config.h>
          #if !(INCLUDE_ONEBODY>=0 && INCLUDE_ERI>=0 && INCLUDE_ERI3>=0 && INCLUDE_ERI2>=0)
          # error \"Libint2 library must be configured with support for 1-body and (2, 3, and 4-center) 2-body integrals\"
          #endif
          int main(int argc, char** argv) {
            return 0;
          }
      "  LIBINT2_IS_CONFIGURED_CORRECTLY)    
    if (NOT LIBINT2_IS_CONFIGURED_CORRECTLY)
      message(FATAL_ERROR "Libint2 library is not configured correctly: need support for 1-body and (2, 3, and 4-center) 2-body integrals")
    endif()
    
    cmake_pop_check_state()
  
    set(LIBINT2_FOUND TRUE)
    set(LIBINT2_LIBRARIES ${LIBINT2_LIBRARY})
    set(LIBINT2_INCLUDE_DIRS ${LIBINT2_INCLUDE_DIR})
    set(LIBINT2_LIBRARY_DIRS ${LIBINT2_LIBRARY_DIR})
    set(LIBINT2_EXTRA_DEFINITIONS -DLIBINT2_DOES_NOT_INLINE_ENGINE)
    mark_as_advanced(LIBINT2_INCLUDE_DIRS LIBINT2_LIBRARY_DIRS LIBINT2_LIBRARIES)
    message(STATUS "Found Libint library: ${LIBINT2_LIBRARIES}")
  endif(LIBINT2_INCLUDE_DIR AND LIBINT2_LIBRARY)


endif(LIBINT2_INCLUDE_DIRS)
