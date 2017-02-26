set(sources
        direct_task_integrals.h
        direct_tile.h
        f12_utility.cpp
        f12_utility.h
        integrals.h
        integral_builder.h
        libint2_engine.cpp
        make_engine.h
        task_integral_kernels.h
        task_integral_kernels.cpp
        screening/screen_base.h
        screening/screen_base.cpp
        screening/schwarz_screen.h
        screening/schwarz_screen.cpp
        task_integrals_common.h
        task_integrals.h
)

add_mpqc_library(lcao_integrals sources sources "MPQClcao_expression;MPQClcao_basis;MPQCmath_groups" "mpqc/chemistry/qc/lcao/integrals")

