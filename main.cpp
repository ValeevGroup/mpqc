#include <iostream>
#include "include/eigen.h"
#include "molecule/cluster_concept.h"
#include "cluster_collapse.h"
#include "molecule/cluster.h"
#include "molecule/molecule.h"
#include "molecule/clustering_functions.h"
#include <vector>
#include "matrix/dmhm.h"
#include <chrono>
#include "molecule/Atom.h"

int main() {
  Atom a({0.0,0.1,0.2},1.0,1.0);
  Atom b({1.0,0.1,0.2},1.0,1.0);
  Atom c({9.0,0.1,0.2},1.0,1.0);
  Atom d({10.0,0.1,0.2},1.0,1.0);
  Cluster ca;
  ca.add_clusterable(c);
  ca.add_clusterable(d);
  ca.guess_center();
  std::cout << ca.center().transpose() << std::endl;

  Molecule mol({a,b,ca});

  auto clusters = mol.cluster_molecule(clustering::kmeans(4),3);

  for(const auto &cluster : clusters){
    std::cout << cluster.center().transpose() << std::endl;
  }
  std::cout << std::endl;

  for(const auto &cluster : clusters){
    for(const auto &atom : collapse_to_atoms(cluster)){
      std::cout << atom.center().transpose() << std::endl;
    }
  }

  for(const auto &cluster : clusters){
    std::cout << cluster.sum_distances_from_center() << std::endl;
  }

  return 0;
}
