set(sources
        ccsd.h
        ccsd.cpp
        ccsd_t.h
        ccsd_t.cpp
        dbccsd.h
        dbccsd.cpp
        diis_ccsd.h
        )

add_library(MPQCcc ${sources})
target_link_libraries(MPQCcc MPQCintegrals MPQCmisc)