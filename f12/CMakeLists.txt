set(sources
    f12_intermediates.h
    f12_utility.cpp
    f12_utility.h
    mp2f12.h
    ccsdf12.h
    f12.cpp
)

add_library(f12 ${sources})
target_link_libraries(f12 ${TA_LIBRARIES} integrals)
