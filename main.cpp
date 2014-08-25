#include <iostream>
#include "include/eigen.h"
#include "molecule/cluster_concept.h"
#include "cluster_collapse.h"
#include "molecule/cluster.h"
#include "molecule/molecule.h"
#include "molecule/clustering_functions.h"
#include <vector>
#include <chrono>
#include "molecule/Atom.h"
#include <Accelerate/Accelerate.h>
#include "include/tbb.h"
#include "include/tiledarray.h"

int main(int argc, char** argv) {
  madness::World &world = madness::initialize(argc, argv);
  world.rank();

  int nthreads = 1;
  std::cout << "input nthreads:";
  std::cin >> nthreads;
  tbb::task_scheduler_init init(nthreads);
  unsigned long N = 5000000;

  tbb::tick_count a0 = tbb::tick_count::now();
  std::vector<Clusterable> atoms;
  atoms.reserve(N);
  for (auto i = 0ul; i < N; ++i) {
    atoms.emplace_back(Atom({0, 0, i}, 1.0, 1.0));
  }
  tbb::tick_count a1 = tbb::tick_count::now();
  double a_alloc = (a1-a0).seconds();
  std::cout << "Atom allocing time = " << a_alloc << std::endl;

  tbb::tick_count m0 = tbb::tick_count::now();
  Molecule mol(std::move(atoms));
  tbb::tick_count m1 = tbb::tick_count::now();
  double m_alloc = (m1-m0).seconds();
  std::cout << "mol allocing time = " << m_alloc << std::endl;

  tbb::tick_count mc0 = tbb::tick_count::now();
  auto clusters = mol.cluster_molecule(clustering::kmeans(10), 60);
  tbb::tick_count mc1 = tbb::tick_count::now();
  double mc_alloc = (mc1-mc0).seconds();
  std::cout << "cluster allocing time = " << mc_alloc << std::endl;

  using iter_t = decltype(clusters.begin());
  double sum = tbb::parallel_reduce(
      tbb::blocked_range<iter_t>(clusters.begin(), clusters.end()), 0.0,
      [](const tbb::blocked_range<iter_t> & r, double d)->double {
        return std::accumulate(r.begin(), r.end(), d,
                               [](double d, const Clusterable &b) {
          return d + b.center().norm();
        });
      },
      std::plus<double>());

  std::cout << "Sum of distances = " << sum << std::endl;

  madness::finalize();
  return 0;
}
