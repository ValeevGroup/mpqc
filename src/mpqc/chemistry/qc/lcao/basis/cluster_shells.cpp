#include "mpqc/chemistry/qc/lcao/basis/cluster_shells.h"
#include "mpqc/chemistry/molecule/cluster.h"

#include <numeric>

namespace mpqc {
namespace lcao {

ClusterShells::ClusterShells() = default;
ClusterShells::ClusterShells(ClusterShells const &) = default;
ClusterShells::ClusterShells(ClusterShells &&) = default;
ClusterShells &ClusterShells::operator=(ClusterShells const &) = default;
ClusterShells &ClusterShells::operator=(ClusterShells &&) = default;

ClusterShells::ClusterShells(std::vector<std::vector<libint2::Shell>> shell,
                             std::shared_ptr<Cluster> c)
    : shells_{std::move(shell)}, cluster_{c} {}

unsigned int ClusterShells::max_am() const { return shells_.size() - 1; }
bool ClusterShells::has_am(unsigned int am) const {
  return shells_.size() > am;
}

unsigned int ClusterShells::nfunctions(unsigned int am) const {
  return std::accumulate(shells_[am].begin(), shells_[am].end(), 0,
                         [&](unsigned int nfuncs, libint2::Shell const &sh) {
                           return nfuncs + sh.size();
                         });
}

unsigned int ClusterShells::flattened_nfunctions() const {
  return std::accumulate(
      shells_.begin(), shells_.end(), 0,
      [&](unsigned int nfuncs, std::vector<libint2::Shell> const &shs) {
        return nfuncs + std::accumulate(
                            shs.begin(), shs.end(), 0,
                            [&](unsigned int nfuncs, libint2::Shell const &sh) {
                              return nfuncs + sh.size();
                            });
      });
}

std::vector<libint2::Shell> const &ClusterShells::shells(
    unsigned int am) const {
  return shells_[am];
}

unsigned int ClusterShells::nshells(unsigned int am) const {
  return shells_[am].size();
}

unsigned int ClusterShells::nshells() const {
  return std::accumulate(
      shells_.begin(), shells_.end(), 0u,
      [](unsigned int n, std::vector<libint2::Shell> const &am_set) {
        return n + am_set.size();
      });
}

std::vector<libint2::Shell> ClusterShells::flattened_shells() const {
  std::vector<libint2::Shell> shells;
  shells.reserve(nshells());

  for (auto const &a : shells_) {
    shells.insert(shells.end(), a.begin(), a.end());
  }

  return shells;
}

Cluster const &ClusterShells::cluster() const { return *cluster_; }

}  // namespace lcao
}  // namespace mpqc
