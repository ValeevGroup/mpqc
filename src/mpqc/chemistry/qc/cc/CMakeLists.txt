set(sources
        ccsd.h
        ccsd_t.h
        dbccsd.h
        diis_ccsd.h
        eom/eom_ccsd.h
        eom/eom_ccsd.cpp
        ccsd.cpp
        ccsd_t.cpp
        )

add_library(MPQCcc ${sources})
target_link_libraries(MPQCcc MPQCintegrals MPQCmisc)