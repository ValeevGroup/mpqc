include(AddLibrary)

# ref library

add_mpqc_object_library(ref
  ref.cc
  reftestx.cc)

# tests

add_executable(reftest EXCLUDE_FROM_ALL reftest.cc $<TARGET_OBJECTS:ref>)
add_test(reftest reftest)
target_link_libraries(
  reftest
  util
)

