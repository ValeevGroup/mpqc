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
#include "include/tbb.h"

int main() {
  tbb::task_scheduler_init init(4);
  unsigned long N = 100000;

  std::vector<Clusterable> atoms;
  atoms.reserve(N);
  for (auto i = 0ul; i < N; ++i) {
    atoms.emplace_back(Atom({ 0, 0, i }, 1.0, 1.0));
  }

  Molecule mol(atoms);

  auto clusters = mol.cluster_molecule(clustering::kmeans(10), 1000);

  using iter_t = decltype(clusters.begin());
  double sum = tbb::parallel_reduce(
      tbb::blocked_range<iter_t>(clusters.begin(), clusters.end()), 0.0,
      [](const tbb::blocked_range<iter_t> & r, double d)->double {
        return std::accumulate(r.begin(), r.end(), d,
                               [](double d, const Cluster &b) {
          return d + b.sum_distances_from_center();
        });
      },
      std::plus<double>());

  std::cout << "Sum of distances = " << sum << std::endl;

  return 0;
}
