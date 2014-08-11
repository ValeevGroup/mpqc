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

  Clusterable ca(a);
  Clusterable cb(b);
  std::cout << ca.center().transpose() << std::endl;
  std::cout << cb.center().transpose() << std::endl;
  using std::swap;
  swap(ca,cb);
  std::cout << ca.center().transpose() << std::endl;
  std::cout << cb.center().transpose() << std::endl;

  Molecule mol({ca,cb});



  return 0;
}
