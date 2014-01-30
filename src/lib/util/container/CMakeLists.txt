include(AddLibrary)

# container library

add_mpqc_object_library(container
  avl.cc
  bitarray.cc
  carray.cc)

# tests

if (MPQC_UNITTEST)
  set_property(SOURCE avltest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(avltest avltest.cc $<TARGET_OBJECTS:container>)
  target_link_libraries(avltest
    util
  )
  add_test(avltest avltest)
endif()

