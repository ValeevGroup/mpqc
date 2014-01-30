include(AddLibrary)

# state library

add_mpqc_object_library(state
  proxy.cc
  state.cc
  state_bin.cc
  state_file.cc
  state_text.cc
  statein.cc
  stateout.cc
  tmplinst.cc
  translate.cc)

# tests

if (MPQC_UNITTEST)
  set_property(SOURCE statetest.cc
    PROPERTY COMPILE_DEFINITIONS SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(statetest statetest.cc $<TARGET_OBJECTS:state>)
  target_link_libraries(statetest
    util
  )
  add_test(statetest statetest)
endif()

