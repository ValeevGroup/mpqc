# optimize library

add_library(optimize OBJECT
  conv.cc
  dfp.cc
  diis.cc
  efc.cc
  function.cc
  gdiis.cc
  kvopt.cc
  lbfgs.f
  mcsearch.cc
  newton.cc
  opt.cc
  powell.cc
  qnewton.cc
  scextrap.cc
  scextrapmat.cc
  steep.cc
  transform.cc
  update.cc)

# tests

set_property(SOURCE opttest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(opttest EXCLUDE_FROM_ALL opttest.cc $<TARGET_OBJECTS:optimize>)
add_test(opttest opttest)

set_property(SOURCE scextest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(scextest EXCLUDE_FROM_ALL scextest.cc $<TARGET_OBJECTS:optimize>)
add_test(scextest scextest)

