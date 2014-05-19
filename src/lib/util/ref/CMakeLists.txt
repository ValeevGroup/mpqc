include(AddLibrary)

# ref library

add_mpqc_object_library(ref
  ref.cc
  reftestx.cc)

# tests

if (MPQC_UNITTEST)
  add_executable(reftest reftest.cc $<TARGET_OBJECTS:ref>)
  target_link_libraries(reftest
    util
  )
  add_test(reftest reftest)
endif()

