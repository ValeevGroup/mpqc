# madness library

add_library(madness OBJECT
  init.cc
)

# tests

#set_property(SOURCE madtest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(madtest EXCLUDE_FROM_ALL madtest.cc $<TARGET_OBJECTS:madness>)
#add_test(madtest madtest)

