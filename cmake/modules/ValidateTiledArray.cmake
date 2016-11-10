include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

macro (validate_tiledarray)

  cmake_push_check_state()

  list(APPEND CMAKE_REQUIRED_LIBRARIES tiledarray)

  ##############################################
  # sanity check: try compiling a simple program
  ##############################################
  CHECK_CXX_SOURCE_COMPILES(
  "
  #include <tiledarray.h>
  int main(int argc, char** argv) {
    TA::World& world = TA::initialize(argc, argv);
    TA::finalize();
    return 0;
  }
  "  TILEDARRAY_COMPILES)

  if (NOT TILEDARRAY_COMPILES)
    message(FATAL_ERROR "TiledArray found, but does not compile correctly.")
  endif()

  ##########################
  # ensure it's fresh enough
  ##########################
  set (TILEDARRAY_OLDEST_REVISION 7087fd4b)
  CHECK_CXX_SOURCE_COMPILES(
  "
  #include <tiledarray.h>
  int main(int argc, char** argv) {
    TA::World& world = TA::initialize(argc, argv);
    
    // test 1
    {
      TA::TArrayD arr;
      auto x = TA::array_to_eigen<TA::Tensor<double>, TA::DensePolicy, Eigen::RowMajor>(arr);
    }
    
    // add more tests here
    
    TA::finalize();
    return 0;
  }
  "  TILEDARRAY_IS_FRESH)

  if (NOT TILEDARRAY_IS_FRESH)
    message(FATAL_ERROR "TiledArray found, but is not fresh enough. Use ${TILEDARRAY_OLDEST_REVISION} or later")
  endif()

  cmake_pop_check_state()
  
endmacro (validate_tiledarray)
