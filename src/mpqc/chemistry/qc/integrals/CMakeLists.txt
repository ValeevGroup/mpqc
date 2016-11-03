set(integral_sources
        screening/screen_base.h
        screening/schwarz_screen.h
        screening/schwarz_screen.cpp
        screening/cached_shell_info.h
        screening/cached_shell_info.cpp
        screening/qqr_screening.h
        screening/qqr_screening.cpp
        screening/qvl_screening.h
        screening/qvl_screening.cpp
        screening/qvl_shell_info.h
        screening/qvl_shell_info.cpp
        ao_factory.h
        ao_factory_base.h
        ao_factory_base.cpp
        direct_ao_factory.h
        lcao_factory.h
        integrals.h
        integral.cpp
        task_integral_kernels.h
        task_integral_kernels.cpp
        libint2_engine.cpp
        periodic_ao_factory.h
)

add_library(MPQCintegrals ${integral_sources})
target_link_libraries(MPQCintegrals ${LIBINT2_LIBRARIES} ${TBB_LIBRARIES} MPQCexpression MPQCbasis MPQCf12utility)

