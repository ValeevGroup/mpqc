include(AddLibrary)

# dft library

add_mpqc_object_library(dft
  am05.cc
  clks.cc
  functional.cc
  hsosks.cc
  integrator.cc
  lebedev.c
  solvent.cc
  tmplinst.cc
  uks.cc
)

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE dfttest.cc PROPERTY COMPILE_DEFINITIONS
              SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(dfttest dfttest.cc $<TARGET_OBJECTS:dft>)
 add_test(dfttest dfttest)
 target_link_libraries( dfttest chemistry )
 
endif()

