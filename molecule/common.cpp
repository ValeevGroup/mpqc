
#include "common.h"
#include "../include/tbb.h"


position_t center_of_mass(const std::vector<Clusterable> cs, double mass){
  tbb::spin_mutex myMutex;
  position_t center = {0, 0, 0};
  tbb::parallel_for(tbb::blocked_range<unsigned long>(0, cs.size()),
                    [&](const tbb::blocked_range<unsigned long> &r) {
    position_t local_sum = {0, 0, 0};
    for (unsigned long z = r.begin(); z != r.end(); ++z) {
      local_sum += cs[z].center() * cs[z].mass();
    }
    tbb::spin_mutex::scoped_lock lock(myMutex);
    center += local_sum;
  });
  center /= mass;
  return center;
}
