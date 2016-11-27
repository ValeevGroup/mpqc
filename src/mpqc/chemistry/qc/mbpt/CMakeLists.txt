set(sources
        dbmp2.h
        dbmp2.cpp
        mp2.h
        mp2.cpp
        )

add_mpqc_library(mbpt sources sources "MPQCscf" "chemistry/qc/mbpt")
