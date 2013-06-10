# solvent library

add_library(solvent OBJECT
  bem.cc
  disprep.cc
)

# tests

#set_property(SOURCE bemtest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(bemtest EXCLUDE_FROM_ALL bemtest.cc $<TARGET_OBJECTS:solvent>)
#add_test(bemtest bemtest)

