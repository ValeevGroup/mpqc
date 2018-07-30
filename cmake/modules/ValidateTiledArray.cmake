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
  set (TILEDARRAY_OLDEST_REVISION 6d8d538769de0c33f6b2bf076b459d54fdf76200)
  # first make sure TiledArray_Eigen and TiledArray_BTAS targets are defined
  if (NOT TARGET TiledArray_Eigen OR NOT TARGET TiledArray_BTAS)
    message (FATAL_ERROR "TiledArray found, but is not fresh enough. Use ${TILEDARRAY_OLDEST_REVISION} or later")
  endif()
  # basic version check
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
    
    // test 2
    {
      auto w = TA::push_default_world(world);
      TA::TArrayD arr;
      auto x = arr(\"i,j\").set_world(world);
    }
    
    // add more tests here
    
    TA::finalize();
    return 0;
  }
  "  TILEDARRAY_IS_FRESH)
  # some version checks will result in failed compilation ... to avoid confusing users
  set(CMAKE_REQUIRED_QUIET_ ${CMAKE_REQUIRED_QUIET})
  set(CMAKE_REQUIRED_QUIET TRUE)
  CHECK_CXX_SOURCE_COMPILES(
          "
  #include <tiledarray.h>
  int main(int argc, char** argv) {
    TA::World& world = TA::initialize(argc, argv);

    // test 1: multiply expression cannot be implicitly converted to a scalar
    {
      TA::TArrayD arr;
      double x = arr(\"i,j\") * arr(\"i,j\");
    }

    // add more tests here

    TA::finalize();
    return 0;
  }
  "  SHOULD_FAIL_IF_TILEDARRAY_IS_FRESH)
  set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_})
  if (SHOULD_FAIL_IF_TILEDARRAY_IS_FRESH)
    message(STATUS "additional checks for TiledArray version failed!")
  endif()
  if (NOT TILEDARRAY_IS_FRESH OR SHOULD_FAIL_IF_TILEDARRAY_IS_FRESH)
    message(FATAL_ERROR "TiledArray found, but is not fresh enough. Use ${TILEDARRAY_OLDEST_REVISION} or later")
  endif()

  ##########################
  # ensure madness::World::get_default is disabled
  ##########################
  CHECK_CXX_SOURCE_COMPILES(
  "
  #include <tiledarray.h>
  #ifndef MADNESS_DISABLE_WORLD_GET_DEFAULT
  # error \"madness::World::get_default is usable\"
  #endif
  int main(int argc, char** argv) {
    return 0;
  }
  "  MADNESS_WORLD_GET_DEFAULT_IS_DISABLED)

  if (NOT MADNESS_WORLD_GET_DEFAULT_IS_DISABLED)
    message(FATAL_ERROR "TiledArray uses external MADNESS, requires MADNESS configured with -DDISABLE_WORLD_GET_DEFAULT=ON CMake option")
  endif()

  cmake_pop_check_state()
  
endmacro (validate_tiledarray)
