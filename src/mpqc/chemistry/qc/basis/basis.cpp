// #include "../molecule/cluster.h"

#include "mpqc/chemistry/molecule/molecule.h"

#include "mpqc/chemistry/qc/basis/basis.h"
#include "mpqc/chemistry/qc/basis/shell_vec_functions.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("Basis", mpqc::basis::Basis);

namespace mpqc {

namespace basis {

Basis::Basis() = default;
Basis::~Basis() = default;
Basis::Basis(Basis const &) = default;
Basis::Basis(Basis &&) = default;

Basis &Basis::operator=(Basis const &) = default;
Basis &Basis::operator=(Basis &&) = default;

Basis::Basis(std::vector<ShellVec> shells) : shells_(std::move(shells)) {}

Basis::Basis(const KeyVal &kv) {
  //  std::cout << "Construct Basis " << std::endl;
  // name of basis
  std::string basis_name = kv.value<std::string>("name");

  // construct basis set
  BasisSet basis_set(basis_name);

  // molecule
  auto mol_ptr = kv.keyval("molecule").class_ptr<mpqc::Molecule>();

  // find world from one level above
  madness::World *world = kv.value<madness::World *>("$:world");

  auto basis = parallel_construct_basis(*world, basis_set, *mol_ptr);

  std::size_t reblock_size = kv.value<int>("reblock",0);
  if(reblock_size != 0){
    basis =  reblock(basis, basis::reblock_basis, reblock_size);
  }

  shells_ = std::move(basis.shells_);
}

int64_t Basis::nfunctions() const {
  int64_t nfuncs = 0;
  for (auto const &shellvec : shells_) {
    for (auto const &sh : shellvec) {
      nfuncs += sh.size();
    }
  }
  return nfuncs;
}

TiledArray::TiledRange1 Basis::create_trange1() const {
  auto blocking = std::vector<int64_t>{0};
  for (auto const &shell_vec : shells_) {
    auto next = blocking.back() + basis::nfunctions(shell_vec);
    blocking.emplace_back(next);
  }

  return TiledArray::TiledRange1(blocking.begin(), blocking.end());
}

int64_t Basis::max_nprim() const {
  int64_t max = 0;
  for (auto const &shell_vec : shells_) {
    const auto current = basis::max_nprim(shell_vec);
    max = std::max(current, max);
  }
  return max;
}

int64_t Basis::max_am() const {
  int64_t max = 0;
  for (auto const &shell_vec : shells_) {
    const auto current = basis::max_am(shell_vec);
    max = std::max(current, max);
  }
  return max;
}

int64_t Basis::nshells() const {
  return std::accumulate(
      shells_.begin(), shells_.end(), int64_t(0),
      [](int64_t x, ShellVec const &a) { return x + int64_t(a.size()); });
}

std::vector<Shell> Basis::flattened_shells() const {
  std::vector<Shell> shells;
  shells.reserve(nshells());

  for (auto const &cluster : cluster_shells()) {
    for (auto const &shell : cluster) {
      shells.push_back(shell);
    }
  }

  return shells;
}

std::vector<ShellVec> const &Basis::cluster_shells() const { return shells_; }

std::ostream &operator<<(std::ostream &os, Basis const &b) {
  unsigned int n = 0;
  for (auto const &shell_vec : b.cluster_shells()) {
    os << "Cluster " << n << "\n";
    ++n;

    for (auto const &shell : shell_vec) {
      os << shell << "\n";
    }
  }

  return os;
}

Basis parallel_construct_basis(madness::World &world, const BasisSet &basis_set,
                               const mpqc::Molecule &mol) {
  Basis basis;

  if (world.rank() == 0) {
    basis = Basis(basis_set.get_cluster_shells(mol));
  }

  world.gop.broadcast_serializable(basis, 0);

  return basis;
}

Basis Basis::join(const Basis &basis) {
  auto self_shells = this->cluster_shells();
  auto other_shells = basis.cluster_shells();
  self_shells.insert(self_shells.end(), other_shells.begin(),
                     other_shells.end());

  return basis::Basis(self_shells);
}
}  // namespace basis
}  // namespace mpqc
