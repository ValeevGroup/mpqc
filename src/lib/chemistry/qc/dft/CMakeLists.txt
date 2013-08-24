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

#set_property(SOURCE dfttest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(dfttest EXCLUDE_FROM_ALL dfttest.cc $<TARGET_OBJECTS:dft>)
#add_test(dfttest dfttest)

