#include <cluster.h>
#include <functional>
#include "../include/tbb.h"
#include "common.h"
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

void Cluster::guess_center() { center_ = center_of_mass(elements_, mass_); }

double Cluster::sum_distances_from_center() const {
    auto reduce_r = [&](double d, const Clusterable &c) {
        return d + std::sqrt(diff_squaredNorm(c.center(), center_));
    };

    using iter_t = decltype(elements_.begin());
    tbb::affinity_partitioner ap;

    return tbb::parallel_reduce(
        tbb::blocked_range<iter_t>(elements_.begin(), elements_.end()), 0.0,
        [&](const tbb::blocked_range<iter_t> &r, double d) {
            return std::accumulate(r.begin(), r.end(), d, reduce_r);
        },
        std::plus<double>(), ap);
}
