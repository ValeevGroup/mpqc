# intv3 library

add_library(intv3 OBJECT
  array.cc
#  atominfo.cc
  bounds.cc
  build2e.cc
  comp1e.cc
  comp2e.cc
  comp2e3c.cc
  init2e.cc
  int1e.cc
  int2e.cc
  intv3.cc
  obintv3.cc
  offsets.cc
  print2e.cc
  shift2e.cc
  storage.cc
  store.cc
  tbintv3.cc
  tformv3.cc
)

# tests

#set_property(SOURCE fjttest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(fjttest EXCLUDE_FROM_ALL fjttest.cc $<TARGET_OBJECTS:intv3>)
#add_test(fjttest fjttest)

#set_property(SOURCE inttest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(inttest EXCLUDE_FROM_ALL inttest.cc $<TARGET_OBJECTS:intv3>)
#add_test(inttest inttest)

