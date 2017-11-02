set(sources
        ccsd.h
        ccsd_hbar.h
        ccsd.cpp
        ccsd_t.h
        ccsd_t.cpp
        eom/eom_ccsd.h
        eom/eom_ccsd_impl.h
        eom/eom_ccsd.cpp
        eom/ea_eom_ccsd.h
        eom/ea_eom_ccsd_impl.h
        eom/ip_eom_ccsd.h
        eom/ip_eom_ccsd_impl.h
        gamma_point_ccsd.h
        gamma_point_ccsd.cpp
        laplace_transform.h
        linkage.h
        solvers.h
        )

add_mpqc_library(lcao_cc sources sources "MPQCcc;MPQClcao_scf;MPQCutil_misc;MPQCmath_quadrature" "mpqc/chemistry/qc/lcao/cc")
