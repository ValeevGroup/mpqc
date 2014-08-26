#include <iostream>
#include "../include/eigen.h"
#include "cluster_concept.h"
#include "cluster_collapse.h"
#include "cluster.h"
#include "molecule.h"
#include "clustering_functions.h"
#include <vector>
#include <chrono>
#include "Atom.h"
#include "../include/tbb.h"

int main(int argc, char** argv) {
  int nthreads = (argc > 1) ? std::stoi(argv[1]) : 4;
  unsigned long N = (argc > 2) ? std::stoi(argv[2]) : 5000;

  tbb::task_scheduler_init init(nthreads);

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

  return 0;
}

