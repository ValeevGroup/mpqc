# ref library

add_library(ref OBJECT
  identity.cc
  ref.cc
  reftestx.cc)

# tests

set_property(SOURCE reftest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(reftest EXCLUDE_FROM_ALL reftest.cc $<TARGET_OBJECTS:ref>)
add_test(reftest reftest)

