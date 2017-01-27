set(sources
    energy.cpp
    energy.h
    findif.cpp
    findif.h
    linkage.h
    property.h
)

add_mpqc_library(properties sources sources "MPQCwfn" "mpqc/chemistry/qc/properties")
