include(AddLibrary)

# keyval library

add_mpqc_object_library(keyval
  ipv2.cc
  ipv2_alloc.cc
  ipv2_cwk.cc
  ipv2_data.cc
  ipv2_error.cc
  ipv2_karray.cc
  ipv2_parse.cc
  ipv2_print.cc
  ipv2_read.cc
  ipv2_scan.cc
  keyval.cc
  keyvalagg.cc
  keyvalass.cc
  keyvalipv2.cc
  keyvalpre.cc
  keyvalstr.cc
  keyvalval.cc)

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE keyvaltest.cc PROPERTY COMPILE_DEFINITIONS
      SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(keyvaltest keyvaltest.cc $<TARGET_OBJECTS:keyval>)
  target_link_libraries(keyvaltest
    util
  )
  add_test(keyvaltest keyvaltest)
endif()

