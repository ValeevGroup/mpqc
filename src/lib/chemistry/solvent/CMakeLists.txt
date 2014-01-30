include(AddLibrary)

# solvent library

add_mpqc_object_library(solvent
  bem.cc
  disprep.cc
)

# tests

#set_property(SOURCE bemtest.cc PROPERTY COMPILE_DEFINITIONS
#    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(bemtest EXCLUDE_FROM_ALL bemtest.cc $<TARGET_OBJECTS:solvent>)
#add_test(bemtest bemtest)

