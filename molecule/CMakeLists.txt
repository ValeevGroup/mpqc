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
