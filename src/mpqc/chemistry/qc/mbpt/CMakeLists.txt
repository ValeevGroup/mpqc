set(sources
        mp2.h
        mbpt.cpp
        )

add_library(MPQCmbpt ${sources})
target_link_libraries(MPQCmbpt ${TA_LIBRARIES} integrals)