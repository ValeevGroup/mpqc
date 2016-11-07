set(sources
        atom.h
        atom.cpp
        atom_based_cluster_concept.h
        atom_based_cluster.h
        atom_based_cluster.cpp
        cluster_concept.h
        cluster.h
        cluster.cpp
        cluster_collapse.h
        clustering_functions.h
        clustering_functions.cpp
        common.h
        common.cpp
        formula.h
        formula.cpp
        molecule.h
        molecule_fwd.h
        molecule.cpp
        unit_cell.h
        unit_cell.cpp
)

add_library(MPQCmolecule ${sources})
target_link_libraries(MPQCmolecule MPQCutil_mad ${TBB_LIBRARIES} ${TA_LIBRARIES})

add_executable(molecule_test EXCLUDE_FROM_ALL examples/molecule_test.cpp)
target_link_libraries(molecule_test MPQCmolecule)
