include(AddLibrary)

# molecule library

add_mpqc_object_library(molecule
  atom.cc
  atominfo.cc
  bend.cc
  cartcoor.cc
  coor.cc
  deriv.cc
  energy.cc
  findisp.cc
  formula.cc
  frag.cc
  imcoor.cc
  linip.cc
  linop.cc
  molden.cc
  molecule.cc
  molfreq.cc
  molrender.cc
  molshape.cc
  molsymm.cc
  out.cc
  redund.cc
  simple.cc
  stors.cc
  stre.cc
  symmcoor.cc
  taylor.cc
  tors.cc
)

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE moltest.cc PROPERTY COMPILE_DEFINITIONS
               SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}" TOPSRCDIR="${CMAKE_SOURCE_DIR}")
  add_executable(moltest moltest.cc $<TARGET_OBJECTS:molecule>)
  add_test(moltest moltest)
  target_link_libraries( moltest chemistry )
endif()

