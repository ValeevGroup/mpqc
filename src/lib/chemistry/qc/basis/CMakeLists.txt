include(AddLibrary)

# basis library

# collect the .cc and .cpp files to be compiled.
SET(ccfiles 
  aotoso.cc
  cart.cc
  cartiter.cc
  dercent.cc
  distsh.cc
  distshpair.cc
  extent.cc
  files.cc
  fjt.cc
  gaussbas.cc
  gaussbaval.cc
  gaussshell.cc
  gaussshval.cc
  gpetite.cc
  integral.cc
  intparams.cc
  inttraits.cc
  lselect.cc
  obint.cc
  operator.cc
  orthog.cc
  petite.cc
  shellrot.cc
  sobasis.cc
  sointegral.cc
  split.cc
  symmint.cc
  tbint.cc
  transform.cc
  uncontract.cc
  union.cc
)

IF(MPQC_NEW_FEATURES)
    set(newccfiles
      tiledbasisset.cpp
    )

    # Compiler Warnings for New Code
    if(FUTURE_WARNING_BUILD)
      foreach(f ${newccfiles})
        set_source_files_properties(f COMPILE_FLAGS "-Wall -Wextra ")
      endforeach()
    endif(FUTURE_WARNING_BUILD)

    LIST(APPEND ccfiles ${newccfiles})
ENDIF(MPQC_NEW_FEATURES)

add_mpqc_object_library(basis ${ccfiles})

IF(MPQC_NEW_FEATURES)
  set_source_files_properties(integral_scaling_test.cpp EXCLUDE_FROM_ALL COMPILE_FLAGS "-Wall -Wextra ")
  set_property(SOURCE integral_scaling_test.cpp PROPERTY COMPILE_DEFINITIONS
        SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(integral_scaling_test integral_scaling_test.cpp $<TARGET_OBJECTS:basis>)
  target_link_libraries(integral_scaling_test chemistry extern_moinfo)
ENDIF(MPQC_NEW_FEATURES)

# tests
IF(MPQC_UNITTEST)
  set_property(SOURCE btest.cc PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(btest btest.cc $<TARGET_OBJECTS:basis>)
  add_test(btest btest)
  
  target_link_libraries( btest chemistry )
  
  set_property(SOURCE tiledbasisset_test.cpp PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(tiledbasisset_test tiledbasisset_test.cpp $<TARGET_OBJECTS:basis>)
  add_test(tiledbasisset_test tiledbasisset_test)
  
  target_link_libraries( tiledbasisset_test chemistry )

  set_property(SOURCE shellorder_test.cpp PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(shellorder_test shellorder_test.cpp $<TARGET_OBJECTS:basis>)
  
  target_link_libraries( shellorder_test chemistry )

  set_property(SOURCE kcluster_test.cpp PROPERTY COMPILE_DEFINITIONS
    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_executable(kcluster_test kcluster_test.cpp $<TARGET_OBJECTS:basis>)
  
  target_link_libraries( kcluster_test chemistry )
  
endif()
