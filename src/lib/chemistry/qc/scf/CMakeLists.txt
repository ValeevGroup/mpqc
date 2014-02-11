include(AddLibrary)

# scf library

# collect the .cc and .cpp files to be compiled.
SET(ccfiles 
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

IF(MPQC_NEW_FEATURES)
    list(APPEND ccfiles tascf.cpp)
ENDIF(MPQC_NEW_FEATURES)

add_mpqc_object_library(scf ${ccfiles})

# tests
if (MPQC_UNITTEST)
  set_property(SOURCE focktest.cc PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(focktest focktest.cc $<TARGET_OBJECTS:basis>)
  add_test(focktest focktest)
  
  target_link_libraries( focktest chemistry )

  set_property(SOURCE scftest.cc PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(scftest scftest.cc $<TARGET_OBJECTS:basis>)
  add_test(scftest scftest)
  
  target_link_libraries( scftest chemistry )

  set_property(SOURCE tascf_test.cpp PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(tascf_test tascf_test.cpp $<TARGET_OBJECTS:basis>)
  add_test(tascf_test tascf_test)
  
  target_link_libraries( tascf_test chemistry )

endif (MPQC_UNITTEST)

