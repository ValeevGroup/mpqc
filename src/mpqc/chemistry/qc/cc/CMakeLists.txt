set(sources
        ccsd.h
        ccsd.cpp
        ccsd_t.h
        ccsd_t.cpp
        dbccsd.h
        dbccsd.cpp
        diis_ccsd.h
        )

add_mpqc_library(cc sources sources "MPQCscf;MPQCutil_misc" "chemistry/qc/cc")
