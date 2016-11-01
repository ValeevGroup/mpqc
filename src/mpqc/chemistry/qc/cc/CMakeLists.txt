set(sources
        ccsd.h
        ccsd_intermediates.h
        ccsd_t.h
        dbccsd.h
        diis_ccsd.h
        eom/eom_ccsd.h
        eom/eom_ccsd.cpp
        integral_generator.h
        lazy_tile.h
        mo_block.h
        cc.cpp
        )

add_library(MPQCcc ${sources})
target_link_libraries(MPQCcc ${TA_LIBRARIES} MPQCintegrals MPQCmisc)