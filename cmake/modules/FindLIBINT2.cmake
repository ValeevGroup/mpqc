
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

if (NOT ${EIGEN3_FOUND})
  find(Eigen3)
endif()

if (LIBINT_INCLUDE_DIRS)
  set (LIBINT_FOUND TRUE)

else (LIBINT_INCLUDE_DIRS)
  set (LIBINT_FOUND FALSE)

  if(LIBINT_INSTALL_DIR)
    set(_LIBINT_INSTALL_DIR ${LIBINT_INSTALL_DIR})
  endif(LIBINT_INSTALL_DIR)

  set(LIBINT_INC_SEARCH_DIR ${_LIBINT_INSTALL_DIR}/include)

  find_path(LIBINT_INCLUDE_DIR NAMES libint2.hpp
    HINTS
    ${LIBINT_INC_SEARCH_DIR})

  find_path(LIBINT_INCLUDE_LIBINT2_DIR NAMES cxxapi.h
    HINTS
    ${LIBINT_INC_SEARCH_DIR}/libint2)

  set(LIBINT_INCLUDE_DIR ${LIBINT_INCLUDE_DIR} ${LIBINT_INCLUDE_LIBINT2_DIR})

  mark_as_advanced(LIBINT_INCLUDE_DIR)

  set(_LIBINT_LIB_NAMES
    "int2"
    )

  set(_LIBINT_LIBRARY_DIR ${_LIBINT_INSTALL_DIR}/lib)
  set(LIBINT_LIBRARY "")

  foreach(_lib ${_LIBINT_LIB_NAMES})
    set(current_lib NOTFOUND)
    find_library(current_lib ${_lib} HINTS ${_LIBINT_LIBRARY_DIR})
    if(NOT current_lib)
      message(FATAL_ERROR "MPQC could not find Libint lib: ${_lib}")
    endif(NOT current_lib)
    LIST(APPEND LIBINT_LIBRARY ${current_lib})
  endforeach()

  set(LIBINT_LIBRARY_DIR ${_LIBINT_LIBRARY_DIR})
  mark_as_advanced(LIBINT_LIBRARY)

  # CODA
  if(LIBINT_INCLUDE_DIR AND LIBINT_LIBRARY)
  
    # validate version, etc. by compiling tests
    cmake_push_check_state()
  
    list(APPEND CMAKE_REQUIRED_INCLUDES ${EIGEN3_INCLUDE_DIR} ${LIBINT_INCLUDE_DIR})
    list(APPEND CMAKE_REQUIRED_LIBRARIES ${LIBINT_LIBRARY})
    
    # sanity check
    CHECK_CXX_SOURCE_COMPILES(
      "
      #include <libint2.hpp>
      int main(int argc, char** argv) {
        libint2::initialize();
        libint2::finalize();
        return 0;
      }
      "  LIBINT_COMPILES)
    if (NOT LIBINT_COMPILES)
      message(FATAL_ERROR "Could not compile LIBINT test program.\nSee CMakeFiles/CMakeError.log for details")
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
            // need libint2::Engine and libint2::TwoBodyEngine
            libint2::Engine e0;
            libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb> e1;
            return 0;
          }
      "  LIBINT_IS_UP_TO_DATE)    
    if (NOT LIBINT_IS_UP_TO_DATE)
      message(FATAL_ERROR "Libint library is too old: 2.2.0 (alpha or later) is required")
    endif()
    
    cmake_pop_check_state()
  
    set(LIBINT_FOUND TRUE)
    set(LIBINT_LIBRARIES ${LIBINT_LIBRARY})
    set(LIBINT_INCLUDE_DIRS ${LIBINT_INCLUDE_DIR})
    set(LIBINT_LIBRARY_DIRS ${LIBINT_LIBRARY_DIR})
    mark_as_advanced(LIBINT_INCLUDE_DIRS LIBINT_LIBRARY_DIRS LIBINT_LIBRARIES)
    message(STATUS "Found Libint libraries: ${LIBINT_LIBRARIES}")
  endif(LIBINT_INCLUDE_DIR AND LIBINT_LIBRARY)


endif(LIBINT_INCLUDE_DIRS)
