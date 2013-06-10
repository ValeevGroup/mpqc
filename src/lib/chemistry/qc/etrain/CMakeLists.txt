# etrain library

add_library(etrain OBJECT
  etrain.cc
)

# tests

#set_property(SOURCE etraintest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(etraintest EXCLUDE_FROM_ALL etraintest.cc $<TARGET_OBJECTS:etrain>)
#add_test(etraintest etraintest)

