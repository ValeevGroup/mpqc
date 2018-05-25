set(cc_sources
        ccsd.h
        ccsd_hbar.h
        ccsd_intermediates.h
        ccsd_r1_r2.h
        ccsd.cpp
        ccsd_t.h
        ccsd_t.cpp
        eom/eom_ccsd.h
        eom/eom_ccsd_impl.h
        eom/eom_preconditioner.h
        eom/eom_ccsd.cpp
        eom/eom_ea_ccsd.h
        eom/eom_ea_ccsd_impl.h
        eom/eom_ip_ccsd.h
        eom/eom_ip_ccsd_impl.h
        gamma_point_ccsd.h
        gamma_point_ccsd.cpp
        laplace_transform.h
        linkage.h
        solvers.h
        )

add_mpqc_library(lcao_cc cc_sources cc_sources "MPQCcc;MPQClcao_scf;MPQCutil_misc;MPQCmath_quadrature" "mpqc/chemistry/qc/lcao/cc")

set(ccsdt_sources
        ccsdt1.h
        ccsdt1.cpp
        ccsdt1b.h
        ccsdt1b.cpp
        ccsdt2.h
        ccsdt2.cpp
        ccsdt3.h
        ccsdt3.cpp
        ccsdt4.h
        ccsdt4.cpp
        ccsdt.h
        ccsdt.cpp
        cc3.h
        cc3.cpp
        ccsdt_linkage.h
        solvers.h
        )

add_mpqc_library(lcao_ccsdt ccsdt_sources ccsdt_sources "MPQCcc;MPQClcao_scf;MPQCutil_misc" "mpqc/chemistry/qc/lcao/cc")
