include(AddLibrary)

# solvent library

add_mpqc_object_library(solvent
  bem.cc
  disprep.cc
)

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE bemtest.cc PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(bemtest bemtest.cc $<TARGET_OBJECTS:basis>)
  add_test(bemtest bemtest)
  
  target_link_libraries(bemtest chemistry)
endif(MPQC_UNITTEST)
