set(sources
        dbmp2.h
        dbmp2.cpp
        denom.h
        linkage.h
        mp2.h
        mp2.cpp
        )

add_mpqc_library(mbpt sources sources "MPQCscf" "mpqc/chemistry/qc/mbpt")
