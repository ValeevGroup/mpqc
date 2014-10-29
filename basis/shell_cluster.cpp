#include "shell_cluster.h"

namespace tcc {
namespace basis {

ShellCluster::ShellCluster() = default;
ShellCluster::ShellCluster(ShellCluster const &) = default;
ShellCluster::ShellCluster(ShellCluster &&) = default;
ShellCluster &ShellCluster::operator=(ShellCluster const &) = default;
ShellCluster &ShellCluster::operator=(ShellCluster &&) = default;

ShellCluster::ShellCluster(std::vector<libint2::Shell> shell,
                           molecule::Cluster c)
    : shells_{std::move(shell)},
      cluster_(std::make_shared<molecule::Cluster>(std::move(c))) {}

std::vector<libint2::Shell> const &ShellCluster::shells() const {
    return shells_;
}

molecule::Cluster const &ShellCluster::cluster() const { return *cluster_; }

} // namespace basis
} // namespace tcc
