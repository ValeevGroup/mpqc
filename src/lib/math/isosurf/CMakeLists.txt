include(AddLibrary)

# isosurf library

add_mpqc_object_library(isosurf
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

if (MPQC_UNITTEST)
  set_property(SOURCE isotest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(isotest isotest.cc $<TARGET_OBJECTS:isosurf>)
  target_link_libraries(isotest
    math
  )
  add_test(isotest isotest)
endif()

