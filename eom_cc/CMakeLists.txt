set(sources 
    eom_ccsd.h
    eom_ccsd.cpp
        )

add_library(libeomcc ${sources})
target_link_libraries(libeomcc ${TA_LIBRARIES})
