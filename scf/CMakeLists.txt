set(sources 
    scf.h
    scf.cpp
    orbital_localization.h
    orbital_localization.cpp
        )

add_library(scf ${sources})
target_link_libraries(scf ${TA_LIBRARIES})
