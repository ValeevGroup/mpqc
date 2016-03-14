set(sources 
    scf.h
    scf.cpp
    orbital_localization.h
    orbital_localization.cpp
    eigen_solve_density_builder.h
    eigen_solve_density_builder.cpp
    purification_density_build.h
    purification_density_build.cpp
        )

add_library(scf ${sources})
target_link_libraries(scf ${TA_LIBRARIES} tensor ta_routines)
