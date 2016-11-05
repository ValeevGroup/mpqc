set(sources
        dbmp2.h
        dbmp2.cpp
        mp2.h
        mp2.cpp
        )

add_library(MPQCmbpt ${sources})
target_link_libraries(MPQCmbpt MPQCwfn ${TA_LIBRARIES} MPQCintegrals)