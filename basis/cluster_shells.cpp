#include "cluster_shells.h"

namespace tcc {
namespace basis {

ClusterShells::ClusterShells() = default;
ClusterShells::ClusterShells(ClusterShells const &) = default;
ClusterShells::ClusterShells(ClusterShells &&) = default;
ClusterShells &ClusterShells::operator=(ClusterShells const &) = default;
ClusterShells &ClusterShells::operator=(ClusterShells &&) = default;

ClusterShells::ClusterShells(std::vector<std::vector<libint2::Shell>> shell,
                             std::shared_ptr<molecule::Cluster> c)
    : shells_{std::move(shell)}, cluster_{c} {}

std::vector<libint2::Shell> const &ClusterShells::shells(amg_mo am) const {
    return shells_[am];
}

molecule::Cluster const &ClusterShells::cluster() const { return *cluster_; }

} // namespace basis
} // namespace tcc
