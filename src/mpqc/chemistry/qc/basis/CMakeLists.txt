set(sources
basis_registry.cpp
basis_registry.h
basis_set_maps.h
basis.cpp
basis.h
cluster_shells.cpp
cluster_shells.h
linkage.h
shell_vec_functions.cpp
shell_vec_functions.h
)

add_mpqc_library(basis sources sources "MPQCmolecule;MPQCmath_ta" "mpqc/chemistry/qc/basis")
target_include_directories(MPQCbasis PUBLIC ${LIBINT2_INCLUDE_DIRS})
target_link_libraries(MPQCbasis PUBLIC ${LIBINT2_LIBRARIES})
target_compile_definitions(MPQCbasis PUBLIC ${LIBINT2_EXTRA_DEFINITIONS})
