set(sources
basis_registry.cpp
basis_registry.h
basis.cpp
basis.h
cluster_shells.cpp
cluster_shells.h
linkage.h
shell_vec_functions.cpp
shell_vec_functions.h
shift_basis.h
shift_basis.cpp
)

add_mpqc_library(lcao_basis sources sources "MPQCmolecule;MPQClcao_expression;MPQCmath_ta" "mpqc/chemistry/qc/lcao/basis")
target_include_directories(MPQClcao_basis PUBLIC ${LIBINT2_INCLUDE_DIRS})
target_link_libraries(MPQClcao_basis PUBLIC ${LIBINT2_LIBRARIES})
target_compile_definitions(MPQClcao_basis PUBLIC ${LIBINT2_EXTRA_DEFINITIONS})
