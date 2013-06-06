# misc library

add_library(misc OBJECT
  bug.cc
#  ccaenv.cc
  compute.cc
  consumableresources.cc
  exenv.cc
  formio.cc
  formiot.cc
  ieee.cc
  regtime.cc
  runnable.cc
  scexception.cc
  tregtime.cc
  units.cc)

# tests

set_property(SOURCE autovectest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(autovectest EXCLUDE_FROM_ALL autovectest.cc $<TARGET_OBJECTS:misc>)
add_test(autovectest autovectest)

set_property(SOURCE bugtest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(bugtest EXCLUDE_FROM_ALL bugtest.cc $<TARGET_OBJECTS:misc>)
add_test(bugtest bugtest)

set_property(SOURCE comptest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(comptest EXCLUDE_FROM_ALL comptest.cc $<TARGET_OBJECTS:misc>)
add_test(comptest comptest)

set_property(SOURCE scinttest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(scinttest EXCLUDE_FROM_ALL scinttest.cc $<TARGET_OBJECTS:misc>)
add_test(scinttest scinttest)

set_property(SOURCE unittest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(unittest EXCLUDE_FROM_ALL unittest.cc $<TARGET_OBJECTS:misc>)
add_test(unittest unittest)

