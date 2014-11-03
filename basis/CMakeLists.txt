set(sources
  basis_set.h
  basis_set.cpp
  atom_basisset.h
  atom_basisset.cpp
  basis_set_maps.h
)

add_library(basis ${sources})

add_executable(libintshell_test EXCLUDE_FROM_ALL libintshell_test.cpp)
target_link_libraries(libintshell_test ${LIBINT_LIBRARIES} basis)
