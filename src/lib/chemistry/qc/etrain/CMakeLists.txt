include(AddLibrary)

# etrain library

add_mpqc_object_library(etrain
  etrain.cc
)

# tests

#set_property(SOURCE etraintest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(etraintest EXCLUDE_FROM_ALL etraintest.cc $<TARGET_OBJECTS:etrain>)
#add_test(etraintest etraintest)

