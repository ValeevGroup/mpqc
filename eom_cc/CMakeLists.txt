set(sources 
    eom_cc.h
    eom_cc.cpp
        )

add_library(libeomcc ${sources})
target_link_libraries(libeomcc ${TA_LIBRARIES})
