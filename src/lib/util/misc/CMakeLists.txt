include(AddLibrary)

# misc library

add_mpqc_object_library(misc
  bug.cc
  compute.cc
  consumableresources.cc
  exception.cc
  exenv.cc
  formio.cc
  ieee.cc
  regtime.cc
  runnable.cc
  scexception.cc
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

set_property(SOURCE formiotest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(formiotest EXCLUDE_FROM_ALL formiotest.cc $<TARGET_OBJECTS:misc>)
add_test(formiotest formiotest)

set_property(SOURCE regtimetest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(regtimetest EXCLUDE_FROM_ALL regtimetest.cc $<TARGET_OBJECTS:misc>)
add_test(regtimetest regtimetest)

set_property(SOURCE scinttest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(scinttest EXCLUDE_FROM_ALL scinttest.cc $<TARGET_OBJECTS:misc>)
add_test(scinttest scinttest)

set_property(SOURCE unittest.cc PROPERTY COMPILE_DEFINITIONS
    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
add_executable(unittest EXCLUDE_FROM_ALL unittest.cc $<TARGET_OBJECTS:misc>)
add_test(unittest unittest)

