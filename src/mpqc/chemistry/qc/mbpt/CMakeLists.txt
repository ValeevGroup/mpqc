set(sources
        dbmp2.h
        linkage.cpp
        mp2.h
        mp2.cpp
        )

add_library(MPQCmbpt ${sources})
target_link_libraries(MPQCmbpt ${TA_LIBRARIES} MPQCintegrals)