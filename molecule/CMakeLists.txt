set(sources
  atom.h
  atom.cpp
  cluster_concept.h
  cluster.h
  cluster.cpp
  cluster_collapse.h
  atom_based_cluster_concept.h
  atom_based_cluster.h
  atom_based_cluster.cpp
  molecule.h
  molecule_fwd.h
  clustering_functions.h
  clustering_functions.cpp
  common.h
  common.cpp
  molecule.cpp
)

add_library(molecule ${sources})
target_link_libraries(molecule ${TBB_LIBRARIES} ${TA_LIBRARIES})

add_executable(molecule_test EXCLUDE_FROM_ALL examples/molecule_test.cpp)
target_link_libraries(molecule_test molecule)
