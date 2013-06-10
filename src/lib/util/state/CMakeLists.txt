# state library

add_library(state OBJECT
  proxy.cc
  state.cc
  state_bin.cc
  state_file.cc
  state_text.cc
  statein.cc
  stateout.cc
  tmplinst.cc
  translate.cc)

# tests

set_property(SOURCE statetest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(statetest EXCLUDE_FROM_ALL statetest.cc $<TARGET_OBJECTS:state>)
add_test(statetest statetest)

