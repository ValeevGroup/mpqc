set(sources
    cc_utility.h
    cc_utility.cpp
    trange1_engine.h
    trange1_engine.cpp
)

add_library(utility ${sources})
target_link_libraries(utility ${TA_LIBRARIES})

