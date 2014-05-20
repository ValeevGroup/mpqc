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

IF(MPQC_NEW_FEATURES AND MPQC_HAS_ELEMENTAL)
    set(newccfiles
      tascf.cpp
      taclscf.cpp
      taclhf.cpp
      cldfgengine.cpp
    )

    # Compiler Warnings for new code
  if(Future_warning_build)
    foreach(f ${newccfiles})
      set_source_files_properties(f COMPILE_FLAGS "-Wall -Wextra ")
    endforeach()
  endif(Future_warning_build)

    list(APPEND ccfiles ${newccfiles})

ENDIF(MPQC_NEW_FEATURES AND MPQC_HAS_ELEMENTAL)

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
  
  
  add_executable(gengine_test cldfgengine_test.cpp $<TARGET_OBJECTS:basis>)
  add_test(gegine_test gengine_test)
  
  target_link_libraries( gengine_test chemistry)

endif (MPQC_UNITTEST)

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/benzene_trimer.kv.in
    ${CMAKE_CURRENT_BINARY_DIR}/benzene_trimer.kv
    @ONLY
  )

if(MPQC_NEW_FEATURES AND MPQC_HAS_ELEMENTAL)
  add_executable(david_hack david_hack_file.cpp $<TARGET_OBJECTS:basis>)
  target_link_libraries(david_hack chemistry)
endif(MPQC_NEW_FEATURES AND MPQC_HAS_ELEMENTAL)

