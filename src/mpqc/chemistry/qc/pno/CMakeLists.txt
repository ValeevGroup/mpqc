set(sources
    ccsd_pno.h
    ccsd_pno.cpp
    ccsd_pno_impl.h
    linkage.h
    )

add_mpqc_library(pno sources sources "MPQCscf" "mpqc/chemistry/qc/pno")
