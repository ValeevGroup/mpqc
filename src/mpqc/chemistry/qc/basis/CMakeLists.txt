set(sources
        basis_set.h
        basis_set.cpp
        basis.h
        basis.cpp
        basis_registry.h
        basis_registry.cpp
        shell_vec_functions.h
        shell_vec_functions.cpp
        cluster_shells.h
        cluster_shells.cpp
        atom_basisset.h
        atom_basisset.cpp
        basis_set_maps.h
)

add_library(MPQCbasis ${sources})
target_link_libraries(MPQCbasis MPQCmolecule utility)

add_executable(libint2_basis_test EXCLUDE_FROM_ALL libint2_basis_test.cpp)
target_link_libraries(libint2_basis_test ${LIBINT2_LIBRARIES} MPQCbasis)
