set(integral_helper_sources
    task_integrals_helper.h
    task_integrals_helper.cpp
    integral_screening_matrices.h
    integral_screening_matrices.cpp
    )
add_library(integrals ${integral_helper_sources})

