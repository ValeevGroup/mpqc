set(sources
  basis_set.h
  basis_set.cpp
  basis.h
  basis.cpp
  shell_vec_functions.h
  shell_vec_functions.cpp
  cluster_shells.h
  cluster_shells.cpp
  atom_basisset.h
  atom_basisset.cpp
  basis_set_maps.h
)

add_library(basis ${sources})
target_link_libraries(basis molecule)

add_executable(libint2_basis_test EXCLUDE_FROM_ALL libint2_basis_test.cpp)
target_link_libraries(libint2_basis_test ${LIBINT_LIBRARIES} basis)
