set(sources
    energy.cpp
    energy.h
    property.cpp
    property.h
)

add_mpqc_library(properties sources sources "MPQCwfn" "mpqc/chemistry/qc/properties")
