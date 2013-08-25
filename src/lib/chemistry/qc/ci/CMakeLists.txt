include(AddLibrary)

# ci library

add_mpqc_object_library(ci
  ci.cc
  compute.cc
)

# tests

#set_property(SOURCE citest.cc PROPERTY COMPILE_DEFINITIONS
#    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(citest EXCLUDE_FROM_ALL citest.cc $<TARGET_OBJECTS:ci>)
#add_test(citest citest)

