set(sources
        ccsd.h
        ccsd_t.h
        dbccsd.h
        diis_ccsd.h
        eom/eom_ccsd.h
        eom/eom_ccsd.cpp
        linkage.cpp
        )

add_library(MPQCcc ${sources})
target_link_libraries(MPQCcc MPQCintegrals MPQCmisc)