set(sources
        linkage.cpp
        linkage.h
        rhf.h
        rhf_impl.h
        rhf.cpp
        orbital_localization.h
        orbital_localization.cpp
        eigen_solve_density_builder.h
        eigen_solve_density_builder_impl.h
        purification_density_build.h
        purification_density_build_impl.h
        util.h
   )

if(TA_DEFAULT_POLICY EQUAL 1)
list(APPEND sources
        zrhf.h
        zrhf.cpp
        pbc/periodic_cond_ortho.h
        pbc/periodic_soad.h
        )
endif()

add_mpqc_library(lcao_scf sources sources "MPQCproperties;MPQClcao_wfn;MPQCmath_ta;MPQCmath_linalg;MPQCmath_clr" "mpqc/chemistry/qc/lcao/scf")
