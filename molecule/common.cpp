
#include "common.h"
#include "../include/tbb.h"
#include "cluster.h"
#include <algorithm>


position_t center_of_mass(const std::vector<Clusterable> cs, double mass) {

  tbb::combinable<position_t> tls_center(position_t({0,0,0}));
  tbb::affinity_partitioner ap;
  tbb::parallel_for(tbb::blocked_range<unsigned long>(0, cs.size()),
                    [&](const tbb::blocked_range<unsigned long> &r) {
    unsigned long z = r.begin();
    unsigned long end = r.end();
    for (; z != end; ++z) {
      tls_center.local() += cs[z].center() * cs[z].mass();
    }
  },ap);

  return tls_center.combine([](
      const position_t &a, const position_t &b) { return a + b; })/mass;
}
