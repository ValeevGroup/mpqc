set(sources
    energy.cpp
    energy.h
    findif.cpp
    findif.h
    linkage.h
    property.h
)

add_mpqc_library(properties sources sources "MPQCwfn;MPQCmath_function;MPQCutil_misc;MPQCutil_c++;tiledarray" "mpqc/chemistry/qc/properties")
