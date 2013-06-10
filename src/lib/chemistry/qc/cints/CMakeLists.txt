# cints library

add_library(cints OBJECT
  cints.cc
  comp_eri.cc
  comp_grt.cc
  edipole.cc
  equadrupole.cc
  eri.cc
  grt.cc
  hcore.cc
  int1e.cc
  int2e.cc
  kinetic.cc
  nuclear.cc
  obintcints.cc
  obosrr.cc
  overlap.cc
  permute2e.cc
  primpairs.cc
  shellpairs.cc
  storage.cc
  tbintcints.cc
  tform.cc
)

# tests

#set_property(SOURCE cintstest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(cintstest EXCLUDE_FROM_ALL cintstest.cc $<TARGET_OBJECTS:cints>)
#add_test(cintstest cintstest)

