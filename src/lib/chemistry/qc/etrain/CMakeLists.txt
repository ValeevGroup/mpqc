include(AddLibrary)

# etrain library

add_mpqc_object_library(etrain
  etrain.cc
)

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE etraintest.cc PROPERTY COMPILE_DEFINITIONS
              SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(etraintest etraintest.cc $<TARGET_OBJECTS:etrain>)
  add_test(etraintest etraintest)
  target_link_libraries( etraintest chemistry )
endif()
  

