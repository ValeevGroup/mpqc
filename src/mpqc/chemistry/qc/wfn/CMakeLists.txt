set(sources
ao_wfn.cpp
ao_wfn.h
lcao_wfn.cpp
lcao_wfn.h
linkage.cpp
linkage.h
spin.cpp
spin.h
trange1_engine.cpp
trange1_engine.h
wfn_forward.h
wfn_world.cpp
wfn_world.h
wfn.cpp
wfn.h
)

add_mpqc_library(wfn sources sources "MPQCintegrals;MPQCutil_keyval;MPQCmolecule" "mpqc/chemistry/qc/wfn")
