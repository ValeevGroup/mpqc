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
  set (TILEDARRAY_OLDEST_REVISION 27622cca)
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

  if (NOT TILEDARRAY_IS_FRESH)
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
