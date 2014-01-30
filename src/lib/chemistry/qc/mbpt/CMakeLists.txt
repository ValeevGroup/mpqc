include(AddLibrary)

# mbpt library

add_mpqc_object_library(mbpt
  bzerofast.c
  cscphf.cc
  csgmat.cc
  csgrad.cc
  csgrad34qb.cc
  csgrade12.cc
  csgrads2pdm.cc
  hsosv1.cc
  hsosv1e1.cc
  hsosv2.cc
  hsosv2lb.cc
  mbpt.cc
  mp2extrap.cc
  util.cc
)

# tests

#set_property(SOURCE mbpttest.cc PROPERTY COMPILE_DEFINITIONS
#    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(mbpttest EXCLUDE_FROM_ALL mbpttest.cc $<TARGET_OBJECTS:mbpt>)
#add_test(mbpttest mbpttest)

