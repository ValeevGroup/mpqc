set(sources
        ccsd.h
        ccsd_intermediates.h
        ccsd_t.h
        diis_ccsd.h
        integral_generator.h
        lazy_tile.h
        mo_block.h
        cc.cpp
        )

add_library(MPQCcc ${sources})
target_link_libraries(MPQCcc ${TA_LIBRARIES} integrals)