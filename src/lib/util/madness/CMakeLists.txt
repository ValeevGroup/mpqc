include(AddLibrary)

# madness library

add_mpqc_object_library(madness
  init.cc world.cc
)

# tests
if (MPQC_UNITTEST)
  add_executable(madtest madtest.cc $<TARGET_OBJECTS:class>)
  target_link_libraries(madtest util)
  add_test(madtest madtest)
endif()

