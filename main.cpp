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
  std::vector<Clusterable> clusterables;

  Atom a(Atom::position_t(0,0,0),16,8);
  clusterables.emplace_back(a);
  for(int i = 0; i < 10; i+=2){
    Cluster temp;
    temp.add_clusterable(Atom(Atom::position_t(i*3,0,i), 1.0, 1.0));
    temp.add_clusterable(Atom(Atom::position_t(i,0,i+1), 1.0, 1.0));
    temp.guess_center();
    clusterables.emplace_back(temp);
  }
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

  std::shuffle(clusterables.begin(), clusterables.end(), std::default_random_engine(seed));

  Molecule mol(clusterables);
  std::vector<Cluster> clusters = mol.cluster_molecule(clustering::kmeans(1), 3);
  for(auto cluster : clusters){
    for(auto a : collapse_to_atoms(cluster)){
      std::cout << a.center().transpose() << std::endl;
    }
    std::cout << std::endl;
  }

  clusterables.emplace_back(mol);
  Molecule mol1(clusterables);
  std::vector<Cluster> clusters1 = mol1.cluster_molecule(clustering::kmeans(1), 3);
  for(auto cluster : clusters1){
    for(auto a : collapse_to_atoms(cluster)){
      std::cout << a.center().transpose() << std::endl;
    }
    std::cout << std::endl;
  }


  return 0;
}
