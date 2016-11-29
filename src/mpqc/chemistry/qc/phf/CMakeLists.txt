set(sources
        linkage.cpp
        phf.h
        phf.cpp
        )

add_mpqc_library(phf sources sources "MPQCmath_ta" "mpqc/chemistry/qc/phf")
