set(sources
        ao_factory.h
        ao_factory_impl.h
        factory_utility.cpp
        factory_utility.h
        factory.cpp
        factory.h
        lcao_factory.h
        linkage.h
        periodic_ao_factory.h
        periodic_lcao_factory.h
        set_oper.h
)

add_mpqc_library(lcao_factory sources sources "MPQClcao_integrals" "mpqc/chemistry/qc/lcao/factory")

