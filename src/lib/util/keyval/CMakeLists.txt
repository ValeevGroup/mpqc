# keyval library

add_library(keyval OBJECT
  ipv2.cc
  ipv2_alloc.cc
  ipv2_cwk.cc
  ipv2_data.cc
  ipv2_error.cc
  ipv2_karray.cc
  ipv2_parse.cc
  ipv2_print.cc
  ipv2_read.cc
  ipv2_scan.cc
  keyval.cc
  keyvalagg.cc
  keyvalass.cc
  keyvalipv2.cc
  keyvalpre.cc
  keyvalstr.cc
  keyvalval.cc)

# tests

set_property(SOURCE keyvaltest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(keyvaltest EXCLUDE_FROM_ALL keyvaltest.cc $<TARGET_OBJECTS:keyval>)
add_test(keyvaltest keyvaltest)

