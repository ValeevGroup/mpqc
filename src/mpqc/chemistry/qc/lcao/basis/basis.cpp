#include "mpqc/chemistry/molecule/common.h"
#include "mpqc/chemistry/molecule/molecule.h"

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/basis/shell_vec_functions.h"
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/util/misc/formio.h"
#include "mpqc/util/misc/exenv.h"

MPQC_CLASS_EXPORT2("Basis", mpqc::lcao::gaussian::AtomicBasis);

namespace mpqc {
namespace lcao {
namespace gaussian {

Basis::Factory::Factory(std::string const &s) : basis_set_name_{s} {}
Basis::Factory::Factory(const KeyVal &kv)
    : basis_set_name_(kv.value<std::string>("name")) {}

std::vector<ShellVec> Basis::Factory::get_cluster_shells(
    Molecule const &mol) const {
  std::vector<ShellVec> cs;
  for (auto const &cluster : mol) {
    const auto libint_atoms =
        ::mpqc::to_libint_atom(collapse_to_atoms(cluster));

    std::streambuf *cout_sbuf = std::cout.rdbuf();  // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    libint2::BasisSet libint_basis(basis_set_name_, libint_atoms);
    std::cout.rdbuf(cout_sbuf);

    // Shells that go with this cluster
    ShellVec cluster_shells;
    cluster_shells.reserve(libint_basis.size());

    for (auto &&shell : libint_basis) {
      cluster_shells.emplace_back(std::move(shell));
    }

    cs.emplace_back(std::move(cluster_shells));
  }

  return cs;
}

ShellVec Basis::Factory::get_flat_shells(Molecule const &mol) const {
  ShellVec cs;
  for (auto const &cluster : mol) {
    const auto libint_atoms = to_libint_atom(collapse_to_atoms(cluster));

    std::streambuf *cout_sbuf = std::cout.rdbuf();  // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    libint2::BasisSet libint_basis(basis_set_name_, libint_atoms);
    std::cout.rdbuf(cout_sbuf);

    for (auto &&shell : libint_basis) {
      cs.emplace_back(std::move(shell));
    }
  }

  return cs;
}

std::ostream& operator<<(std::ostream &os, Basis::Factory const &f) {
  os << indent << "Basis::Factory:\n" << incindent;
  os << indent << "name = " << f.name() << std::endl;
  os << decindent;
  return os;
}

////////////////////

Basis::Basis() = default;
Basis::~Basis() { }
Basis::Basis(Basis const &) = default;
Basis::Basis(Basis &&) = default;

Basis &Basis::operator=(Basis const &) = default;
Basis &Basis::operator=(Basis &&) = default;

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
    auto next = blocking.back() + ::mpqc::lcao::gaussian::nfunctions(shell_vec);
    blocking.emplace_back(next);
  }

  return TiledArray::TiledRange1(blocking.begin(), blocking.end());
}

int64_t Basis::max_nprim() const {
  int64_t max = 0;
  for (auto const &shell_vec : shells_) {
    const auto current = ::mpqc::lcao::gaussian::max_nprim(shell_vec);
    max = std::max(current, max);
  }
  return max;
}

int64_t Basis::max_am() const {
  int64_t max = 0;
  for (auto const &shell_vec : shells_) {
    const auto current = ::mpqc::lcao::gaussian::max_am(shell_vec);
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

std::ostream& operator<<(std::ostream &os, Basis const &b) {
  os << indent << "Basis:\n" << incindent;
  size_t clidx = 0;
  size_t shidx = 0;
  for (auto const &shell_vec : b.cluster_shells()) {
    os << indent << "Cluster " << clidx << ":\n" << incindent;

    for (auto const &shell : shell_vec) {
      os << indent << "Shell " << shidx << ":\n" << incindent;
      os << indent << shell << "\n" << decindent;
      ++shidx;
    }

    os << decindent;

    ++clidx;
  }

  os << decindent;
  return os;
}

Eigen::RowVectorXi sub_basis_map(const Basis &basis, const Basis &sub_basis) {
  auto shells = basis.flattened_shells();
  auto sub_shells = sub_basis.flattened_shells();

  std::size_t n_functions = basis.nfunctions();
  std::size_t n_sub_functions = sub_basis.nfunctions();
  TA_ASSERT(n_functions >= n_sub_functions);

  std::size_t n_shells = basis.nshells();
  std::size_t n_sub_shells = sub_basis.nshells();

  Eigen::RowVectorXi result = Eigen::RowVectorXi::Zero(n_functions);

  std::size_t sub_shell_lowbound = 0;
  std::size_t sub_shell_highbound = 0;
  for (std::size_t i = 0; i < n_sub_shells; i++) {
    auto &sub_shell = sub_shells[i];
    sub_shell_lowbound = sub_shell_highbound;
    sub_shell_highbound += sub_shell.size();

    std::size_t shell_lowbound = 0;
    std::size_t shell_highbound = 0;

    bool find = false;

    // locate the position in shells
    for (std::size_t j = 0; j < n_shells; j++) {
      auto &shell = shells[j];
      shell_lowbound = shell_highbound;
      shell_highbound += shell.size();

      if (sub_shell == shell) {
        find = true;

        // fill in the value
        const std::size_t gap = sub_shell_lowbound - shell_lowbound;
        for (std::size_t k = shell_lowbound; k < shell_highbound; k++) {
          result[k] = gap + k + 1;
        }

        break;
      }
    }

    if (find == false) {
      throw std::runtime_error("\n Not a Sub Basis! \n");
    }
  }

  return result;
}

Basis merge(const Basis &basis1, const Basis &basis2) {
  auto shells1 = basis1.cluster_shells();
  auto shells2 = basis2.cluster_shells();
  shells1.insert(shells1.end(), shells2.begin(), shells2.end());

  return Basis(shells1);
}

Basis parallel_make_basis(madness::World &world, const Basis::Factory &factory,
                          const mpqc::Molecule &mol) {
  Basis basis;
  if (world.rank() == 0) {
    auto shells = factory.get_cluster_shells(mol);
    basis = std::move(shells);
  }

  world.gop.broadcast_serializable(basis, 0);
  return basis;
}

AtomicBasis::AtomicBasis(const KeyVal &kv)
    : factory_(std::make_shared<Factory>(kv)),
      molecule_(kv.class_ptr<mpqc::Molecule>("atoms")) {
  if (!molecule_) {  // use old keyword "molecule"
    molecule_ = kv.class_ptr<mpqc::Molecule>("molecule");
  }
  if (!molecule_)
    throw InputError("AtomicBasis did not receive atoms", __FILE__, __LINE__,
                     "atoms");

  auto *world = kv.value<madness::World *>("$:world");
  static_cast<Basis&>(*this) = parallel_make_basis(*world, *factory_, *molecule_);

  std::size_t reblock_size = kv.value<int>("reblock", 0);
  if (reblock_size != 0) {
    static_cast<Basis&>(*this) = reblock(*this, reblock_basis, reblock_size);
  }

  compute_shell_to_atom();

  auto rebuilder = [this]() { this->rebuild_shells(); };
  Observer::register_message(molecule_.get(), rebuilder);
}

namespace detail {
bool equal(const std::array<double, 3> &arr3,
                const Eigen::Vector3d &vec3) {
  return arr3[0] == vec3(0) && arr3[1] == vec3(1) && arr3[2] == vec3(2);
}
}  // namespace detail

void AtomicBasis::compute_shell_to_atom() {
  const auto ncluster = shells_.size();
  shell_to_atom_.resize(ncluster);

  auto atoms = molecule_->atoms();
  const auto natoms = atoms.size();
  size_t current = 0;
  for (size_t cluster = 0; cluster != ncluster; ++cluster) {
    for (const auto &shell : shells_[cluster]) {
      while (!detail::equal(shell.O,atoms[current].center()) && current != natoms) {
        ++current;
      }
      if (current == natoms)
        throw ProgrammingError("invalid shells_ in AtomicBasis", __FILE__, __LINE__);
      shell_to_atom_[cluster].push_back(current);
    }
  }
  assert(current == natoms-1);
}

void AtomicBasis::rebuild_shells() {
  auto atoms = molecule_->atoms();
  for(size_t c=0; c != shells_.size(); ++c) {
    for(size_t s=0; s!=shells_[c].size(); ++s) {
      auto& shell = shells_[c][s];
      auto atom = shell_to_atom_[c][s];
      const auto& O_new = atoms[atom].center();
      for(int xyz=0; xyz!=3; ++xyz)
        shell.O[xyz] = O_new(xyz);
    }
  }
}

std::shared_ptr<const Molecule>
AtomicBasis::molecule() const { return molecule_; }

std::shared_ptr<const Basis::Factory>
AtomicBasis::factory() const { return factory_; }

std::ostream& operator<<(std::ostream &os, AtomicBasis const &b) {
  os << indent << "AtomicBasis:\n" << incindent;
  os << *b.factory();
  os << *b.molecule();
  os << static_cast<const Basis&>(b);
  os << decindent;
  return os;
}

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
