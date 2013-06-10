# container library

add_library(container OBJECT
  avl.cc
  bitarray.cc
  carray.cc)

# tests

set_property(SOURCE avltest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(avltest EXCLUDE_FROM_ALL avltest.cc $<TARGET_OBJECTS:container>)
add_test(avltest avltest)

