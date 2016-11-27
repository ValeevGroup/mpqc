set(sources
        cabs_singles.h
        f12_intermediates.h
        mp2f12.cpp
        mp2f12.h
        gf2f12.h
        gf2f12.cpp
        ccsd_f12.h
        ccsd_f12.cpp
        ccsd_t_f12.h
        ccsd_t_f12.cpp
        db_f12_intermediates.h
        dbccsd_f12.h
        dbccsd_f12.cpp
        dbmp2f12.h
        dbmp2f12.cpp
        )

add_mpqc_library(f12 sources sources "MPQCscf;MPQCmbpt;MPQCcc" "chemistry/qc/f12")

