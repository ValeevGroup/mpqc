include(AddLibrary)

# psi library
IF(HAVE_PSI3)
    SET(ccfiles
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
    add_mpqc_object_library(psi ${ccfiles})

    # tests
    if (MPQC_UNITTEST)
      set_property(SOURCE psitest.cc PROPERTY COMPILE_DEFINITIONS
                   SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
      add_executable(psitest EXCLUDE_FROM_ALL psitest.cc $<TARGET_OBJECTS:psi>)
      add_test(psitest psitest)
      target_link_libraries(psitest chemistry)
    endif (MPQC_UNITTEST)
    
ENDIF(HAVE_PSI3)
