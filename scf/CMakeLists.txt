set(sources 
    scf.h
    scf.cpp
        )

add_library(scf ${sources})
target_link_libraries(scf ${TA_LIBRARIES})
