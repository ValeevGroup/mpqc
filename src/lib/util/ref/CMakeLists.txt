# ref library

add_library(ref OBJECT
  ref.cc
  reftestx.cc)

# tests

add_executable(reftest EXCLUDE_FROM_ALL reftest.cc $<TARGET_OBJECTS:ref>)
add_test(reftest reftest)
target_link_libraries(
  reftest
  util
)

