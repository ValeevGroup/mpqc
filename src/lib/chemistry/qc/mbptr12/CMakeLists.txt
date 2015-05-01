include(AddLibrary)

# mbptr12 library
IF(MPQC_R12)
    SET(ccfiles
        FxFgen.cc 
        Xgen.cc 
        approxApp.cc 
        approxB.cc 
        approxC.cc 
        approxC_GenRefansatz2.cc 
        approxC_ansatz1.cc 
        approxCp.cc 
        bc_contrib.cc 
        compute_1rdm_X_F12.cc 
        compute_VXB_a.cc 
        compute_VXB_a_vbsneqobs.cc 
        compute_VXB_c_ansatz1.cc 
        compute_VX_GenRefansatz2.cc 
        compute_amps.cc 
        compute_energy_a.cc 
        coulomb.cc 
        creator.cc 
        dkh.cc 
        ebc_contribs.cc 
        exchange.cc
        extern_pt2r12.cc
        fixedcoefficient.cc 
        fock.cc 
        gbc_contribs.cc 
        intermeds_g12.cc 
        mbptr12.cc
        molcas_pt2r12.cc
        mp2_pair_energies.cc 
        mp2r12_energy.cc 
        mp2r12_energy_util.cc 
        pt2r12.cc 
        pt2r12_spinorbital.cc 
        r12_amps.cc 
        r12int_eval.cc 
        r12technology.cc 
        r12wfnworld.cc 
        ri_basis.cc 
        singles_casscf.cc
        singles_emp2.cc 
        sr_r12intermediates_util.cc
        twobody_intermediates_ta.cc 
        twobody_intermeds.cc 
        twobodygrid.cc 
        twoparticlecontraction.cc
    )

    add_mpqc_object_library(mbptr12 ${ccfiles})

ENDIF(MPQC_R12)
# tests

#set_property(SOURCE mbptr12test.cc PROPERTY COMPILE_DEFINITIONS
#    SRCDIR="${CMAKE_CURRENT_SOURCE_DIR}")
#add_executable(mbptr12test EXCLUDE_FROM_ALL mbptr12test.cc $<TARGET_OBJECTS:mbptr12>)
#add_test(mbptr12test mbptr12test)

