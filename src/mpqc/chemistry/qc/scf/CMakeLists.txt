set(sources
        linkage.cpp
        linkage.h
        rhf.h
        rhf.cpp
        orbital_localization.h
        orbital_localization.cpp
        eigen_solve_density_builder.h
        eigen_solve_density_builder.cpp
        purification_density_build.h
        purification_density_build.cpp
        util.h
        pbc/prhf.h
        pbc/prhf.cpp
        )

add_mpqc_library(scf sources sources "MPQCwfn;MPQCmath_ta;MPQCmath_linalg;MPQCmath_clr" "mpqc/chemistry/qc/scf")
