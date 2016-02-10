set(sources
        scf.h
        traditional_df_fock_builder.h
        traditional_df_fock_builder.cpp
        )

add_library(scf ${sources})
target_link_libraries(scf ${TA_LIBRARIES})
