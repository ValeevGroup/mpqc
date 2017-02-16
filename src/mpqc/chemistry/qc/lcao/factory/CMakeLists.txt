set(sources
        ao_factory.h
        factory_utility.h
        factory_utility.cpp
        factory.h
        factory.cpp
        lcao_factory.h
        periodic_ao_factory.cpp
        periodic_ao_factory.h
        periodic_lcao_factory.cpp
        periodic_lcao_factory.h
        set_oper.h
)

add_mpqc_library(lcao_factory sources sources "MPQClcao_integrals" "mpqc/chemistry/qc/lcao/factory")

