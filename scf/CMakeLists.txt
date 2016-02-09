set(sources
        scf.cpp
        scf.h
        )

add_library(scf ${sources})
target_link_libraries(scf ${TA_LIBRARIES})