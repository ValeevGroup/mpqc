include(AddLibrary)

# scf library

add_mpqc_object_library(scf
  clhf.cc
  clscf.cc
  effh.cc
  fbclhf.cc
  hsoshf.cc
  hsosscf.cc
  osshf.cc
  ossscf.cc
  scf.cc
  scfgradient.cc
  scfops.cc
  scfvector.cc
  svd.cc
  tchf.cc
  tcscf.cc
  tmplinst.cc
  uhf.cc
  uscf.cc
)

# tests

#set_property(SOURCE focktest.cc PROPERTY COMPILE_DEFINITIONS
#    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(focktest EXCLUDE_FROM_ALL focktest.cc $<TARGET_OBJECTS:scf>)
#add_test(focktest focktest)

#set_property(SOURCE scftest.cc PROPERTY COMPILE_DEFINITIONS
#    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(scftest EXCLUDE_FROM_ALL scftest.cc $<TARGET_OBJECTS:scf>)
#add_test(scftest scftest)

