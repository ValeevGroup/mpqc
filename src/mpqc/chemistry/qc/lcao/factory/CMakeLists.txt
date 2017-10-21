set(sources
        ao_factory.h
        ao_factory_impl.h
        factory_utility.cpp
        factory_utility.h
        factory.cpp
        factory.h
        lcao_factory.h
        linkage.cpp
        linkage.h
        periodic_ao_factory.cpp
        periodic_ao_factory.h
        periodic_lcao_factory.h
        set_oper.h
        wfn_world.cpp
        wfn_world.h
)

add_mpqc_library(lcao_factory sources sources "MPQClcao_integrals" "mpqc/chemistry/qc/lcao/factory")

