include(AddLibrary)

# render library

add_mpqc_object_library(render
  algebra3.cc
  animate.cc
  appearance.cc
  color.cc
  material.cc
  object.cc
  oogl.cc
  polygons.cc
  polylines.cc
  polysphere.cc
  render.cc
  sphere.cc
  tempinst.cc
  transform.cc)

# tests

if (MPQC_UNITTEST)
  set_property(SOURCE rentest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(rentest rentest.cc $<TARGET_OBJECTS:render>)
  target_link_libraries(rentest
    util
  )
  add_test(rentest rentest)
endif()

