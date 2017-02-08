set(sources
        ../factory/ao_factory.h
        ../factory/factory_utility.h
        ../factory/factory_utility.cpp
        direct_task_integrals.h
        direct_tile.h
        f12_utility.cpp
        f12_utility.h
        integrals.h
        ../factory/factory.cpp
        integral_builder.h
        ../factory/lcao_factory.h
        libint2_engine.cpp
        make_engine.h
        task_integral_kernels.h
        task_integral_kernels.cpp
        screening/screen_base.h
        screening/schwarz_screen.h
        screening/schwarz_screen.cpp
        task_integrals_common.h
        task_integrals.h
)

add_mpqc_library(lcao_integrals sources sources "MPQClcao_expression;MPQClcao_basis" "mpqc/chemistry/qc/lcao/integrals")

