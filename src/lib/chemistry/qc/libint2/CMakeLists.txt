include(AddLibrary)

# libint2 library

IF(HAVE_LIBINT2)
    set(ccfiles
      bounds.cc
      #comp_eri.cc
      comp_g12.cc
      comp_g12dkh.cc
      comp_g12nc.cc
      core_ints_engine.cc
      #comp_grt.cc
      edipole.cc
      efield.cc
      efield_grad.cc
      equadrupole.cc
      #eri.cc
      g12.cc
      g12dkh.cc
      g12nc.cc
      #grt.cc
      gto.cc
      hcore.cc
      int1e.cc
      int2e.cc
      kinetic.cc
      libint2.cc
      nuclear.cc
      obintlibint2.cc
      obosrr.cc
      overlap.cc
      p4.cc
      permute2e.cc
      primpairs.cc
      shellpairs.cc
      tbintlibint2.cc
      tform.cc
    )

    add_mpqc_object_library(libint2 ${ccfiles})

    # tests
    if (MPQC_UNITTEST)
        set_property(SOURCE libint2test.cc PROPERTY COMPILE_DEFINITIONS
                     SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")

        add_executable(libint2test libint2test.cc )
        add_test(libint2test libint2test)
        target_link_libraries(
          libint2test
          chemistry
          ${LIBINT2_LIBRARY}
        )
    endif(MPQC_UNITTEST)

ENDIF(HAVE_LIBINT2)

