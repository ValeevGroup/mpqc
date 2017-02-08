set(sources
        ao_factory.h
        factory_utility.h
        factory_utility.cpp
        factory.cpp
        lcao_factory.h
)

add_mpqc_library(lcao_factory sources sources "MPQClcao_integral" "mpqc/chemistry/qc/lcao/factory")

