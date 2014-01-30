include(AddLibrary)

# madness library

add_mpqc_object_library(madness
  init.cc
)

# tests

#set_property(SOURCE madtest.cc PROPERTY COMPILE_DEFINITIONS
#    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(madtest EXCLUDE_FROM_ALL madtest.cc $<TARGET_OBJECTS:madness>)
#add_test(madtest madtest)

