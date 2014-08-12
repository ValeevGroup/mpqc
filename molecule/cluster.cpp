#include <cluster.h>
#include <functional>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/spin_mutex.h>
Cluster::Cluster(Cluster &&c) noexcept : elements_(std::move(c.elements_)),
                                         charge_(std::move(c.charge_)),
                                         mass_(std::move(c.mass_)) {
  c.center_.swap(center_);
}

Cluster &Cluster::operator=(Cluster &&c) noexcept {
  elements_ = std::move(c.elements_);
  charge_ = std::move(c.charge_);
  mass_ = std::move(c.mass_);
  c.center_.swap(center_);
  return *this;
}

void Cluster::guess_center() {
  tbb::spin_mutex myMutex;
  center_ = { 0, 0, 0 };
  tbb::parallel_for(tbb::blocked_range<unsigned long>(0, elements_.size()),
                    [&](const tbb::blocked_range<unsigned long> &r) {
    position_t local_sum = { 0, 0, 0 };
    for (unsigned long z = r.begin(); z != r.end(); ++z) {
      local_sum += elements_[z].center() * elements_[z].mass();
    }
    tbb::spin_mutex::scoped_lock lock(myMutex);
    center_ += local_sum;
  });
  center_ /= mass_;
}

double Cluster::sum_distances_from_center() const {
  auto reduce_r = [&](double d, const Clusterable &c) {
    return d + (c.center() - center_).norm();
  };

  using iter_t = decltype(elements_.begin());

  return tbb::parallel_reduce(
      tbb::blocked_range<iter_t>(elements_.begin(), elements_.end()), 0.0,
      [&](const tbb::blocked_range<iter_t> &r, double d) {
        return std::accumulate(r.begin(), r.end(), d, reduce_r);
      },
      std::plus<double>());
}
