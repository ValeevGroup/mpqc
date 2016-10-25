set(sources
        f12_utility.cpp
        f12_utility.h
)

add_library(MPQCf12utility ${sources})
target_link_libraries(MPQCf12utility ${TA_LIBRARIES} )

set(sources
        cabs_singles.h
        f12_intermediates.h
        f12_utility.cpp
        f12_utility.h
        mp2f12.cpp
        mp2f12.h
        ccsdf12.h
        db_f12_intermediates.h
        dbccsdf12.h
        dbmp2f12.h
        dbmp2f12.cpp
        f12.cpp
)

add_library(MPQCf12 ${sources})

target_link_libraries(MPQCf12 MPQCintegrals)

