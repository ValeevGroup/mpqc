# psi library

add_library(psi OBJECT
  psicc.cc
  psicc_pt2r12.cc
  psici.cc
  psiexenv.cc
  psifile11.cc
  psiinput.cc
  psiqtorder.cc
  psirdm.cc
  psiref.cc
  psiwfn.cc
)

# tests

#set_property(SOURCE psitest.cc PROPERTY COMPILE_DEFINITIONS
#    -DSRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(psitest EXCLUDE_FROM_ALL psitest.cc $<TARGET_OBJECTS:psi>)
#add_test(psitest psitest)

