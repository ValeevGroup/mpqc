set(sources
        linkage.cpp
        linkage.h
        rhf.h
        rhf_impl.h
        rhf.cpp
        orbital_localization.h
        orbital_localization.cpp
        eigen_solve_density_builder.h
        eigen_solve_density_builder.cpp
        purification_density_build.h
        purification_density_build.cpp
        util.h
        zrhf.h
        zrhf.cpp
        pbc/periodic_cond_ortho.h
        pbc/periodic_soad.h
   )

add_mpqc_library(scf sources sources "MPQCwfn;MPQCmath_ta;MPQCmath_linalg;MPQCmath_clr" "mpqc/chemistry/qc/scf")
