# render library

add_library(render OBJECT
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

set_property(SOURCE rentest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(rentest EXCLUDE_FROM_ALL rentest.cc $<TARGET_OBJECTS:render>)
add_test(rentest rentest)

