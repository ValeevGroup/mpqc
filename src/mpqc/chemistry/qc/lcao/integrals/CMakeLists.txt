set(sources
        density_fitting/cadf_coeffs.h
        density_fitting/cadf_coeffs.cpp
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
        screening/cached_shell_info.cpp
        screening/cached_shell_info.h
        screening/qqr_screening.cpp
        screening/qqr_screening.h
        screening/qvl_screening.cpp
        screening/qvl_screening.h
        screening/qvl_shell_info.cpp
        screening/qvl_shell_info.h
        task_integrals_common.h
        task_integrals.h
)

add_mpqc_library(lcao_integrals sources sources "MPQClcao_expression;MPQClcao_basis;MPQCmath_groups;MPQCutil_misc" "mpqc/chemistry/qc/lcao/integrals")

