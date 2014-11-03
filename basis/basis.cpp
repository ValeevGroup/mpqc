#include "basis.h"
#include "../molecule/cluster.h"
#include "cluster_shells.h"
#include "basis_set.h"

namespace tcc {
namespace basis {

Basis::Basis() = default;
Basis::~Basis() = default;
Basis::Basis(Basis const &) = default;
Basis::Basis(Basis &&) = default;

Basis &Basis::operator=(Basis const &) = default;
Basis &Basis::operator=(Basis &&) = default;

Basis::Basis(std::vector<ClusterShells> cs) : cluster_shells_(std::move(cs)) { }

} // namespace basis
} // namespace tcc
