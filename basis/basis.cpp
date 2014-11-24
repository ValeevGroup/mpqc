#include "basis.h"
#include "../molecule/cluster.h"
#include "cluster_shells.h"
#include "basis_set.h"
#include "../include/libint.h"
#include "../include/tiledarray.h"

#include <vector>

namespace tcc {
namespace basis {

Basis::Basis() = default;
Basis::~Basis() = default;
Basis::Basis(Basis const &) = default;
Basis::Basis(Basis &&) = default;

Basis &Basis::operator=(Basis const &) = default;
Basis &Basis::operator=(Basis &&) = default;

Basis::Basis(std::vector<ClusterShells> cs) : cluster_shells_(std::move(cs)) {}

std::vector<unsigned int> Basis::am_blocking_generator() const {
    std::vector<unsigned int> blocking = {0};
    unsigned int max_am
        = std::max_element(cluster_shells_.begin(), cluster_shells_.end(),
                           [&](ClusterShells const &a, ClusterShells const &b) {
                               return a.max_am() < b.max_am();
                           })->max_am();

    for (auto am = 0u; am < max_am + 1; ++am) {
        for (auto const &cluster : cluster_shells_) {
            if (cluster.has_am(am)) {
                auto next = blocking.back() + cluster.nfunctions(am);
                blocking.emplace_back(next);
            }
        }
    }

    return blocking;
};

std::vector<unsigned int> Basis::flattened_blocking_generator() const {
    std::vector<unsigned int> blocking = {0};

    for (auto const &cluster : cluster_shells_) {
        auto next = blocking.back() + cluster.flattened_nfunctions();
        blocking.emplace_back(next);
    }

    return blocking;
}

TiledArray::TiledRange1 Basis::create_trange1() const {
    auto blocking = am_blocking_generator();
    return TiledArray::TiledRange1{blocking.begin(), blocking.end()};
}

TiledArray::TiledRange1 Basis::create_flattend_trange1() const {
    auto blocking = flattened_blocking_generator();
    return TiledArray::TiledRange1{blocking.begin(), blocking.end()};
}

std::vector<ClusterShells> const &Basis::cluster_shells() const {
    return cluster_shells_;
}

std::ostream &operator<<(std::ostream &os, Basis const &b) {
    unsigned int n = 0;
    for (auto const &cluster : b.cluster_shells()) {
        os << "Cluster " << n << "\n";
        ++n;

        for (auto const &shell : cluster.flattened_shells()) {
            os << shell << "\n";
        }
    }

    return os;
}

} // namespace basis
} // namespace tcc
