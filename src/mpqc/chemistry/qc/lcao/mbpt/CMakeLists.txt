set(sources
        dbmp2.h
        dbmp2.cpp
        dbmp2_impl.h
        denom.h
        linkage.h
        mp2.h
        mp2.cpp
        mp2_impl.h
        gamma_point_mp2.h
        gamma_point_mp2.cpp
        )

add_mpqc_library(lcao_mbpt sources sources "MPQClcao_scf" "mpqc/chemistry/qc/lcao/mbpt")
