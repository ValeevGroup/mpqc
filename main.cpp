#include <iostream>
#include "include/eigen.h"
#include "molecule/cluster_concept.h"
#include "cluster_collapse.h"
#include "molecule/cluster.h"
#include "molecule/molecule.h"
#include "molecule/clustering_functions.h"
#include <vector>
#include "molecule/Atom.h"

int main() {
  std::vector<Clusterable> clusterables;

  for(int i = 0; i < 10; i+=2){
    Cluster temp;
    temp.add_clusterable(Atom(Atom::position_t(0,0,i), 1.0, 1.0));
    temp.add_clusterable(Atom(Atom::position_t(0,0,i+1), 1.0, 1.0));
    temp.guess_center();
    clusterables.emplace_back(temp);
  }

  Molecule mol(clusterables);
  for(int i = 1; i <= clusterables.size(); ++i){
    std::cout << "Number of clusters is " << i << std::endl;
    for(int j = 1; j < 100; ++j){
      std::vector<Cluster> clusters = mol.cluster_molecule(clustering::kmeans(j),i);
      if(min > clustering::sum_cluster_distances(clusters)) {
        for(auto cluster : clusters){
          for(auto a : collapse_to_atoms(cluster)){
            std::cout << a.center().transpose() << std::endl;
          }
          std::cout << std::endl;
        }
        min = clustering::sum_cluster_distances(clusters);
        std::cout << "new min = " << min << std::endl;
      }
    }
    std::cout << std::endl;
  }

  return 0;
}
