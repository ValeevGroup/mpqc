set(sources
        ccsd.h
        ccsd.cpp
        ccsd_t.h
        ccsd_t.cpp
        dbccsd.h
        dbccsd.cpp
        diis_ccsd.h
        gamma_point_ccsd.h
        gamma_point_ccsd.cpp
        linkage.h
        )

add_mpqc_library(cc sources sources "MPQCscf;MPQCmbpt;MPQCutil_misc" "mpqc/chemistry/qc/cc")
