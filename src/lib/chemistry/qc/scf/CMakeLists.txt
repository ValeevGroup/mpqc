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
    list(APPEND ccfiles taclscf.cpp)
    list(APPEND ccfiles taclhf.cpp)
    list(APPEND ccfiles cldfgfactory.cpp)
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

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/tascf_test.kv.in
    ${CMAKE_CURRENT_BINARY_DIR}/tascf_test.kv
    @ONLY
  )
  add_executable(tascf_test tascf_test.cpp $<TARGET_OBJECTS:basis>)
  add_test(tascf_test tascf_test)
  
  target_link_libraries( tascf_test chemistry )
  
  
  add_executable(taclscf_test taclscf_test.cpp $<TARGET_OBJECTS:basis>)
  add_test(taclscf_test taclscf_test)
  
  target_link_libraries( taclscf_test chemistry )

  add_executable(taclhf_test taclhf_test.cpp $<TARGET_OBJECTS:basis>)
  add_test(taclhf_test taclhf_test)
  
  target_link_libraries( taclhf_test chemistry )
  
  
  add_executable(gfactory_test cldfgfactory_test.cpp $<TARGET_OBJECTS:basis>)
  add_test(gfactory_test gfactory_test)
  
  target_link_libraries( gfactory_test chemistry)

endif (MPQC_UNITTEST)

