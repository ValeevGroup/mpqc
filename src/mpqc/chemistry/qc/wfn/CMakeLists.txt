set(sources
    wfn.h
    wfn.cpp
    ao_wfn.h
    ao_wfn.cpp
    lcao_wfn.h
    lcao_wfn.cpp
    wfn_world.h
    wfn_world.cpp
)

add_library(MPQCwfn ${sources})
target_link_libraries(MPQCwfn ${TA_LIBRARIES} MPQCkeyval MPQCmolecule)
