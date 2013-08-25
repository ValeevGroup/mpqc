include(AddLibrary)

# symmetry library

add_mpqc_object_library(symmetry
  chartab.cc
  corrtab.cc
  ico.cc
  irrep.cc
  maketab.cc
  pointgrp.cc
  rep.cc
  symop.cc
  tetra.cc)

# tests
set_property(SOURCE pointgrouptest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(pointgrouptest EXCLUDE_FROM_ALL pointgrouptest.cc $<TARGET_OBJECTS:symmetry>)
add_test(pointgrouptest pointgrouptest)

