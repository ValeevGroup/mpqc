include(AddLibrary)

# intv3 library

set(ccfiles
  array.cc
#  atominfo.cc
  bounds.cc
  build2e.cc
  comp1e.cc
  comp2e.cc
  comp2e3c.cc
  init2e.cc
  int1e.cc
  int2e.cc
  intv3.cc
  obintv3.cc
  offsets.cc
  print2e.cc
  shift2e.cc
  storage.cc
  store.cc
  tbintv3.cc
  tformv3.cc
)

add_mpqc_object_library(intv3 ${ccfiles})

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE fjttest.cc PROPERTY COMPILE_DEFINITIONS
               SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(fjttest fjttest.cc $<TARGET_OBJECTS:intv3>)
  add_test(fjttest fjttest)
  target_link_libraries( fjttest chemistry )

  set_property(SOURCE inttest.cc PROPERTY COMPILE_DEFINITIONS
               SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(inttest inttest.cc $<TARGET_OBJECTS:intv3>)
  add_test(inttest inttest)
  target_link_libraries( inttest chemistry )
endif()


