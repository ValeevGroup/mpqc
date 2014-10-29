set(sources
  atom.h
  cluster_concept.h
  cluster.h
  cluster_collapse.h
  molecule.h
  molecule_fwd.h
  clustering_functions.h
  clustering_functions.cpp
  cluster.cpp
  common.h
  common.cpp
  molecule.cpp
)

add_library(molecule ${sources})
target_link_libraries(molecule ${TBB_LIBRARIES})

add_executable(molecule_test EXCLUDE_FROM_ALL examples/molecule_test.cpp)
target_link_libraries(molecule_test molecule)
