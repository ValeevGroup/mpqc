set(sources
ao_wfn.cpp
ao_wfn.h
lcao_wfn.cpp
lcao_wfn.h
linkage.cpp
linkage.h
spin.cpp
spin.h
wfn_forward.h
wfn_world.cpp
wfn_world.h
wfn.cpp
wfn.h
)

add_mpqc_library(lcao_wfn sources sources "MPQClcao_factory;MPQCutil_keyval;MPQCmolecule;tiledarray" "mpqc/chemistry/qc/lcao/wfn")
