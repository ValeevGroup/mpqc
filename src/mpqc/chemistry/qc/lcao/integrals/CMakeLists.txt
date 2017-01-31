set(sources
        ao_factory.h
        ao_factory_base.h
        ao_factory_base.cpp
        direct_ao_factory.h
        direct_task_integrals.h
        direct_tile.h
        f12_utility.cpp
        f12_utility.h
        integrals.h
        integral.cpp
        integral_builder.h
        lcao_factory.h
        libint2_engine.cpp
        make_engine.h
        periodic_ao_factory.h
        periodic_ao_factory.cpp
        periodic_lcao_factory.h
        task_integral_kernels.h
        task_integral_kernels.cpp
        task_integrals_common.h
        task_integrals.h
        screening/cached_shell_info.h
        screening/cached_shell_info.cpp
        screening/qqr_screening.h
        screening/qqr_screening.cpp
        screening/qvl_screening.h
        screening/qvl_screening.cpp
        screening/qvl_shell_info.h
        screening/qvl_shell_info.cpp
        screening/screen_base.h
        screening/schwarz_screen.h
        screening/schwarz_screen.cpp
)

add_mpqc_library(lcao_integrals sources sources "MPQClcao_expression;MPQClcao_basis" "mpqc/chemistry/qc/lcao/integrals")

