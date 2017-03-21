set(sources
    cis.cpp
    cis.h
    linkage.h
    )

add_mpqc_library(lcao_ci sources sources "MPQClcao_scf;MPQCutil_misc" "mpqc/chemistry/qc/lcao/ci")