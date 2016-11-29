set(sources
    energy.h
    energy.cpp
    propertybase.h
)

add_mpqc_library(properties sources sources "MPQCwfn" "mpqc/chemistry/qc/properties")
