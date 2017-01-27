set(sources
        ao_factory.h
        ao_factory_base.h
        ao_factory_base.cpp
        direct_ao_factory.h
        f12_utility.cpp
        f12_utility.h
        integrals.h
        integral.cpp
        lcao_factory.h
        libint2_engine.cpp
        task_integral_kernels.h
        task_integral_kernels.cpp
        screening/screen_base.h
        screening/schwarz_screen.h
        screening/schwarz_screen.cpp
        # screening/cached_shell_info.h
        # screening/cached_shell_info.cpp
        # screening/qqr_screening.h
        # screening/qqr_screening.cpp
        # screening/qvl_screening.h
        # screening/qvl_screening.cpp
        # screening/qvl_shell_info.h
        # screening/qvl_shell_info.cpp
        periodic_ao_factory.h
        periodic_ao_factory.cpp
)

add_mpqc_library(integrals sources sources "MPQCexpression;MPQCbasis;boost" "mpqc/chemistry/qc/integrals")

