set(sources
    wfn.cpp
    wfn.h
)

add_mpqc_library(wfn sources sources "MPQCmolecule;MPQCutil_misc;MPQCutil_keyval" "mpqc/chemistry/qc/wfn")
