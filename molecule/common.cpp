
#include "common.h"
#include "../include/tbb.h"
#include "cluster.h"
#include <algorithm>


position_t center_of_mass(const std::vector<Clusterable> cs, double mass) {

    // tbb::combinable<position_t> tls_center(position_t({0, 0, 0}));
    // tbb::parallel_for(0ul, cs.size(), [&](unsigned long z) {
    //                                    tls_center.local() +=
    //                                        cs[z].center() * cs[z].mass();
    //                                  },
    //                  );
    return tbb::parallel_reduce(
               tbb::blocked_range<unsigned long>(0, cs.size()),
               position_t({0, 0, 0}),
               [&](const tbb::blocked_range<unsigned long> & r, position_t p)
                   ->position_t {
                   auto i = r.begin();
                   const auto end = r.end();
                   for (; i != end; ++i) {
                       p += cs[i].center() * cs[i].mass();
                   }
                   return p;
               },
               std::plus<position_t>()) / mass;

    // return tls_center.combine([](const position_t &a,
    //                              const position_t &b) { return a + b; }) /
    //        mass;
}
