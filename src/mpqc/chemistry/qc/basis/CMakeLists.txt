set(sources
atom_basisset.cpp
atom_basisset.h
basis_fwd.h
basis_registry.cpp
basis_registry.h
basis_set_maps.h
basis_set.cpp
basis_set.h
basis.cpp
basis.h
cluster_shells.cpp
cluster_shells.h
shell_vec_functions.cpp
shell_vec_functions.h
)

add_mpqc_library(basis sources sources "MPQCmolecule;MPQCmath_ta" "chemistry/qc/basis")
