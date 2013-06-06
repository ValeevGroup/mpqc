# isosurf library

add_library(isosurf OBJECT
  edge.cc
  implicit.c
  isosurf.cc
  shape.cc
  surf.cc
  surfor.cc
  surfse.cc
  surfst.cc
  tmplinst.cc
  triangle.cc
  tricoef.cc
  vertex.cc
  volume.cc
  volume_i.cc)

# tests

set_property(SOURCE isotest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(isotest EXCLUDE_FROM_ALL isotest.cc $<TARGET_OBJECTS:isosurf>)
add_test(isotest isotest)

