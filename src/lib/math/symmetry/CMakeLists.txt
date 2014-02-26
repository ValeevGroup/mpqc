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

if (MPQC_UNITTEST)
  set_property(SOURCE pointgrouptest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(pointgrouptest pointgrouptest.cc $<TARGET_OBJECTS:symmetry>)
  target_link_libraries(pointgrouptest
    math
  )
  add_test(pointgrouptest pointgrouptest)
endif()

