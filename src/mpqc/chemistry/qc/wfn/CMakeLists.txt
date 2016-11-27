set(sources
    wfn.h
    wfn.cpp
    ao_wfn.h
    ao_wfn.cpp
    lcao_wfn.h
    lcao_wfn.cpp
    linkage.cpp
    trange1_engine.h
    trange1_engine.cpp
    wfn_world.h
    wfn_world.cpp
)

add_mpqc_library(wfn sources sources "MPQCintegrals;MPQCutil_keyval;MPQCmolecule" "chemistry/qc/wfn")
