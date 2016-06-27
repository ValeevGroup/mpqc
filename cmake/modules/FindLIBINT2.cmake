
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

if (NOT ${EIGEN3_FOUND})
  find(Eigen3)
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

  find_path(LIBINT2_INCLUDE_LIBINT2_DIR NAMES cxxapi.h
    HINTS
    ${LIBINT2_INC_SEARCH_DIR}/libint2)

  set(LIBINT2_INCLUDE_DIR ${LIBINT2_INCLUDE_DIR} ${LIBINT2_INCLUDE_LIBINT2_DIR})

  mark_as_advanced(LIBINT2_INCLUDE_DIR)

  set(_LIBINT2_LIB_NAMES
    "int2"
    )

  set(_LIBINT2_LIBRARY_DIR ${_LIBINT2_INSTALL_DIR}/lib)
  set(LIBINT2_LIBRARY "")

  foreach(_lib ${_LIBINT2_LIB_NAMES})
    set(current_lib NOTFOUND)
    find_library(current_lib ${_lib} HINTS ${_LIBINT2_LIBRARY_DIR})
    if(NOT current_lib)
      message(FATAL_ERROR "MPQC could not find Libint2 lib: ${_lib}")
    endif(NOT current_lib)
    LIST(APPEND LIBINT2_LIBRARY ${current_lib})
  endforeach()

  set(LIBINT2_LIBRARY_DIR ${_LIBINT2_LIBRARY_DIR})
  mark_as_advanced(LIBINT2_LIBRARY)

  # CODA
  if(LIBINT2_INCLUDE_DIR AND LIBINT2_LIBRARY)
  
    # validate version, etc. by compiling tests
    cmake_push_check_state()
  
    list(APPEND CMAKE_REQUIRED_INCLUDES ${EIGEN3_INCLUDE_DIR} ${LIBINT2_INCLUDE_DIR})
    list(APPEND CMAKE_REQUIRED_LIBRARIES ${LIBINT2_LIBRARY})
    
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
          #if !(LIBINT_MAJOR_VERSION==2 && LIBINT_MINOR_VERSION==2 && LIBINT_MICRO_VERSION>=0)
          # error \"Libint2 library is too old\"
          #endif
          int main(int argc, char** argv) {
            // need libint2::Engine
            libint2::Engine e0;
            return 0;
          }
      "  LIBINT2_IS_UP_TO_DATE)    
    if (NOT LIBINT2_IS_UP_TO_DATE)
      message(FATAL_ERROR "Libint2 library is too old: 2.2.0 (beta) is required")
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
