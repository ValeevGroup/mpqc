set(sources
        ccsd.h
        ccsd.cpp
        ccsd_t.h
        ccsd_t.cpp
        dbccsd.h
        dbccsd.cpp
        gamma_point_ccsd.h
        gamma_point_ccsd.cpp
        linkage.h
        )

add_mpqc_library(lcao_cc sources sources "MPQCcc;MPQClcao_scf;MPQCutil_misc" "mpqc/chemistry/qc/lcao/cc")
